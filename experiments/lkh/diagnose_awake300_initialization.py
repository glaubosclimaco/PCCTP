"""Passive complete-geometry initialization diagnostics; never launches search."""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from experiments.lkh.raw_route_audit import load, check


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == '__main__':
    auditpath = HERE/'ruin_total300_awake_train_audit.json'
    summarypath = HERE/'ruin_total300_awake_train_summary.json'
    envpath = HERE/'ruin_total300_awake_environment_audit.json'
    audit, summary, environment = [json.loads(p.read_text()) for p in (auditpath, summarypath, envpath)]
    assert audit['status'] == environment['status'] == 'PASS'
    assert audit['runs'] == summary['runs']
    snapshot = HERE/'awake300_initialization_snapshots'/f"runs{audit['runs']}_{digest(auditpath)[:12]}"
    snapshot.mkdir(parents=True, exist_ok=True)
    snapshot_hashes = {}
    for source in (auditpath, summarypath, envpath):
        target = snapshot/source.name
        content = source.read_bytes()
        if target.exists(): assert target.read_bytes() == content
        else:
            with target.open('xb') as handle: handle.write(content)
        snapshot_hashes[str(target.relative_to(ROOT))] = digest(target)
    accepted = {(r['instance'], r['seed'], r['mode']): r for r in audit['rows']}
    valid_blocks = {(r['instance'], r['seed']) for r in environment['rows']}
    inputs, rows, hashes, effects = {}, [], {}, []
    for geometry in summary['rows']:
        name = geometry['instance']; data = load(name); inputs[name] = data['input_hash']
        costs = {mode: [] for mode in ('mdm', 'hybrid10', 'hybrid60')}
        for seed in (70, 71, 72):
            assert (name, seed) in valid_blocks
            for mode in costs:
                path = HERE/'ruin_total300_awake_train'/f'{name}.seed{seed}.{mode}.json'
                assert digest(path) == accepted[(name, seed, mode)]['result_sha256']
                r = json.loads(path.read_text()); check(data, r['route'], r['cost'])
                assert r['cost'] == accepted[(name, seed, mode)]['cost']
                costs[mode].append(Fraction(str(r['cost']))); hashes[path.name] = digest(path)
                if mode == 'mdm': continue
                check(data, r['pre_route'], r['pre_cost'])
                native, post = r['native'], r['post']
                assert post['initial_cost'] == r['pre_cost']
                counters = {k: post[k] for k in ('lkh_calls', 'lkh_completed', 'lkh_timeouts', 'lkh_failed')}
                assert counters['lkh_calls'] == sum(counters[k] for k in ('lkh_completed', 'lkh_timeouts', 'lkh_failed'))
                assert (native['cost'] is not None) == native['feasible_before_deadline']
                missing = not native['feasible_before_deadline']
                if missing:
                    assert r['pre_cost'] == r['fallback']['cost']
                    assert r['pre_route'] == r['fallback']['route']
                rows.append(dict(instance=name, seed=seed, mode=mode,
                    no_timely_native=missing, native_cost=native['cost'], initial_cost=r['pre_cost'],
                    final_cost=r['cost'], initial_support=len(r['pre_route']), final_support=len(r['route']),
                    reduction_after_initialization_percent=float(100*(Fraction(str(r['pre_cost']))-Fraction(str(r['cost'])))/Fraction(str(r['pre_cost']))),
                    counters=counters))
        means = {m: sum(v)/3 for m, v in costs.items()}
        gains = {m: 100*(means['mdm']-means[m])/means['mdm'] for m in ('hybrid10', 'hybrid60')}
        for mode, gain in gains.items(): assert abs(float(gain)-geometry['gains'][mode]) < 1e-9
        effects.append(dict(instance=name, gains={m: float(v) for m, v in gains.items()},
            no_timely_native={m: sum(r['no_timely_native'] for r in rows if r['instance']==name and r['mode']==m) for m in gains}))
    result = dict(status='PASS', scope='Passive diagnostics on all currently complete TRAIN geometries; '
                  'no reserved inputs, optimization, new gate, causal claim or independent confirmation.',
                  batch_runs=audit['runs'], batch_complete=audit['batch_complete'], completed_geometries=len(effects),
                  upstream_snapshot_sha256=snapshot_hashes,
                  rows=rows, geometry_effects=effects, input_sha256=inputs, result_sha256=hashes,
                  source_sha256={str(p.relative_to(ROOT)): digest(p) for p in
                                 (Path(__file__), HERE/'raw_route_audit.py', auditpath, summarypath, envpath)},
                  limitations=['Exit code zero does not prove a timely validated improving tour.',
                    'Counter buckets can include the final call after the outer event deadline.',
                    'Other status combines nonzero exit and exhausted setup allowance; no per-call logs separate them.',
                    'Poor initialization and ordering statuses are associated observations, not isolated causal effects.'])
    (HERE/'awake300_initialization_diagnostics.json').write_text(json.dumps(result, indent=2)+'\n')
    lines = ['# Diagnóstico passivo de inicialização no lote300s', '',
             f"Snapshot: {audit['runs']}/81 execuções auditadas; {len(effects)} geometrias completas. Lote completo: {audit['batch_complete']}.", '',
             'Nenhuma alteração de política ou lançamento de otimizador. Todas as geometrias completas são mantidas.', '',
             '| Instância | Ganho H10 (%) | Ganho H60 (%) | Sem solução MDM no prazo H10/H60 |',
             '|---|---:|---:|---:|']
    for e in effects:
        lines.append(f"| {e['instance']} | {e['gains']['hybrid10']:.6f} | {e['gains']['hybrid60']:.6f} | {e['no_timely_native']['hybrid10']}/3; {e['no_timely_native']['hybrid60']}/3 |")
    lines += ['', '## Detalhe de rd400', '',
              '| Semente | Método | MDM tempestivo | Inicial | Final | Suporte inicial | Chamadas/exit0/timeout/other |',
              '|---|---|---|---:|---:|---:|---|']
    for r in rows:
        if not r['instance'].startswith('g_rd400_'): continue
        counts='/'.join(str(r['counters'][k]) for k in ('lkh_calls','lkh_completed','lkh_timeouts','lkh_failed'))
        lines.append(f"| {r['seed']} | {r['mode']} | {not r['no_timely_native']} | {r['initial_cost']:g} | {r['final_cost']:g} | {r['initial_support']} | {counts} |")
    lines += ['', '## Interpretação limitada', '',
              'A inicialização de10s falha em fornecer um incumbente MDM tempestivo em duas sementes de rd400. '
              'A reparação então começa do fallback. Nessas duas fases, o contador exit0 é zero. '
              'Mesmo H60/seed72 tem zero exit0 e ainda reduz o custo: reparo guloso pode produzir melhorias '
              'sem reordenação LKH bem-sucedida. Isso não identifica uma causa única nem prova benefício do cache.', '',
              'A perda de H10 em rd400 é preservada; não excluir a geometria, trocar sementes ou retunar a fase. '
              'O protocolo de300s permanece o mesmo. Comparação MDM70 completa continua necessária.', '']
    counterexample = next((r for r in rows if r['instance'].startswith('g_pr226_') and r['seed']==70
                          and r['mode']=='hybrid10' and r['no_timely_native']), None)
    if counterexample:
        lines += [f"Contraponto: pr226/seed70/H10 também parte do fallback, mas cai de "
                  f"{counterexample['initial_cost']:g} para {counterexample['final_cost']:g}, com "
                  f"{counterexample['counters']['lkh_completed']} registros exit0. A ausência de MDM "
                  'tempestivo, isoladamente, não explica todas as diferenças de qualidade.', '']
    lines += ['- '+text for text in result['limitations']]
    (ROOT/'research_notes/awake300_initialization_diagnostics.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status='PASS', geometries=len(effects), hybrid_records=len(rows), native_optimizer_calls=0)))
