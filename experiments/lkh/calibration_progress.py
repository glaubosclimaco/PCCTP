"""Passive nine-slot progress report; no optimizer or authoritative audit refresh."""
import hashlib
import datetime
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from experiments.lkh.raw_route_audit import load, check


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def observe():
    folder = HERE/'published70_calibration'; protocol_path = folder/'protocol.json'
    assert protocol_path.exists(), 'Calibration has not started'
    protocol = json.loads(protocol_path.read_text())
    assert protocol['seed'] == 80 and protocol['iterations'] == 70 and len(protocol['cases']) == 9
    assert protocol['positive_budget_disabled'] and protocol['safety_timeout_seconds'] == 43200
    roles = json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())['roles']
    assert all(roles[c['geometry']] == 'train' for c in protocol['cases'])
    failure_path = folder/'environment_failure.json'
    failure = json.loads(failure_path.read_text()) if failure_path.exists() else None
    audit_path = HERE/'published70_calibration_audit.json'
    audit = json.loads(audit_path.read_text()) if audit_path.exists() else {}
    verified = {r['instance']: r for r in audit.get('rows', [])} if audit.get('status') == 'PASS' else {}
    rows = []
    for case in protocol['cases']:
        name = case['instance']; target = folder/f'{name}.seed80.json'
        row = dict(instance=name, status='NOT_STARTED', observed_iterations=None,
                   observed_cost=None, reference_duration_seconds=None, independently_audited=False)
        if target.exists():
            content = target.read_bytes(); result = json.loads(content)
            assert result['instance'] == name and result['seed'] == 80
            assert result['status'] in ('COMPLETED70', 'CENSORED_OR_FAILED')
            assert result['native']['budget_seconds'] == 0
            row.update(status=result['status'], observed_iterations=result['observed_iterations'],
                       observed_cost=result['observed_cost'], result_sha256=hashlib.sha256(content).hexdigest())
            if result['status'] == 'COMPLETED70':
                check(load(name), result['route'], result['cost'])
                duration = result['native']['process_wall_seconds']
                assert math.isfinite(duration) and duration > 0
                row['reference_duration_seconds'] = duration
            else: assert result['cost'] is None and result['route'] == []
            row['independently_audited'] = (name in verified and
                    verified[name]['result_sha256'] == row['result_sha256'])
        else:
            checkpoint = folder/f'{name}.seed80.native.raw.checkpoint'
            if checkpoint.exists():
                try:
                    content = checkpoint.read_bytes(); payload = json.loads(content)
                except (OSError, ValueError): row['status'] = 'CHECKPOINT_UNREADABLE_DURING_UPDATE'
                else:
                    assert payload['seed'] == 80 and payload['budget_seconds'] == 0
                    assert type(payload['iterations']) is int and 0 <= payload['iterations'] <= 70
                    check(load(name), payload['route'], payload['cost'])
                    row.update(status='CHECKPOINT_ONLY', observed_iterations=payload['iterations'],
                               observed_cost=payload['cost'], mining_calls=payload['mining_calls'],
                               checkpoint_elapsed_seconds=payload['seconds'],
                               checkpoint_completed_flag=payload['completed'],
                               checkpoint_sha256=hashlib.sha256(content).hexdigest())
            elif (folder/f'{name}.seed80.native.log').exists(): row['status'] = 'RUN_STARTED_NO_CHECKPOINT'
        if failure and failure['instance'] == name:
            row['status'] = 'ENVIRONMENT_FAILURE'; row['environment_failure_sha256'] = digest(failure_path)
        rows.append(row)
    completed = sum(r['status'] == 'COMPLETED70' for r in rows)
    audited = sum(r['independently_audited'] for r in rows)
    full = bool(completed == audited == 9 and audit.get('batch_complete') and audit.get('all_references_completed'))
    projection = 30*sum(r['reference_duration_seconds'] for r in rows)/3600 if full else None
    return dict(status='COMPLETE_AUDITED' if full else 'PROGRESS_ONLY',
        observed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), scheduled_references=9,
        completed70_by_runner=completed, records_covered_by_independent_audit=audited, rows=rows,
        rough_ten_seed_three_method_process_hours=projection,
        scope='Passive calibration progress, not the authoritative final audit or a performance comparison. '
              'Checkpoint costs/iterations are provisional, not70-iteration endpoints. '
              'No linear extrapolation from iteration count. Any full-calibration projection is planning only: '
              '30times the nine reference durations, excluding audit overhead and duration variability.',
        source_sha256={str(p.relative_to(ROOT)):digest(p) for p in (Path(__file__),protocol_path,
                       HERE/'raw_route_audit.py')},
        authoritative_audit_sha256=digest(audit_path) if audit_path.exists() else None)


if __name__ == '__main__':
    report = observe()
    target = HERE/'published70_progress.json'; temporary = target.with_suffix('.tmp')
    temporary.write_text(json.dumps(report,indent=2)+'\n'); temporary.replace(target)
    lines = ['# Progresso passivo da calibração MDM70', '', 'Observado em UTC: '+report['observed_utc'], '',
             f"Estado: {report['status']}. Referências previstas:9; concluídas pelo runner:{report['completed70_by_runner']}; cobertas pelo auditor independente:{report['records_covered_by_independent_audit']}.", '',
             '| Instância | Estado | Iterações observadas | Custo observado (provisório se checkpoint) |',
             '|---|---|---:|---:|']
    for row in report['rows']:
        lines.append(f"| {row['instance']} | {row['status']} | {row['observed_iterations']} | {row['observed_cost']} |")
    lines += ['', 'Custos de checkpoint não são resultados MDM70 completos. Não estimar término por extrapolação linear das iterações. '
              'Este relatório não inicia busca, não altera a fila e não substitui a auditoria final. Nove slots são preservados, incluindo falhas e censuras.', '']
    (ROOT/'research_notes/published70_progress.md').write_text('\n'.join(lines))
    print(json.dumps({k:report[k] for k in ('status','scheduled_references','completed70_by_runner',
                                           'records_covered_by_independent_audit')}))
