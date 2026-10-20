"""Reuse unchanged original observers with explicit fresh-batch output paths."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    from audit_awake_environment import audit
    environment=audit(HERE/'ruin_total300_awake_train')
    (HERE/'ruin_total300_awake_environment_audit.json').write_text(json.dumps(environment,indent=2)+'\n')
    if environment['batch_complete']:assert environment['eligible_complete_environment']
    sources=['audit_ruin_total300_train.py','summarize_ruin_total300_train.py','anytime_curve300.py',
        'profile_mining_exposure.py','report_decision300.py','diagnose_ruin300_regimes.py']
    import sys
    for name in sources:
        path=HERE/name
        source=path.read_text().replace('ruin_total300_train','ruin_total300_awake_train')
        source=source.replace('DECISION300_REPORT.md','DECISION300_AWAKE_REPORT.md')
        source=source.replace('decision300_report_provenance.json','decision300_awake_report_provenance.json')
        source=source.replace('ruin300_regimes.md','ruin300_awake_regimes.md')
        sys.argv=[str(path)]+(['ruin_total300_awake_train'] if name=='profile_mining_exposure.py' else [])
        exec(compile(source,str(path),'exec'),dict(__name__='__main__',__file__=str(path)))
    report_path=HERE.parents[1]/'research_notes/DECISION300_AWAKE_REPORT.md'
    report_path.write_text(f"> Auditoria ambiental independente: {environment['validated_blocks']}/27 blocos validados; "
        f"elegibilidade ambiental do lote completo: {environment['eligible_complete_environment']}. "
        "Execuções no bloco ainda aberto aguardam seu diagnóstico ambiental.\n\n"+report_path.read_text())
    (HERE/'awake300_observer_provenance.json').write_text(json.dumps(dict(
        transformation='literal batch/output-name substitution only; original frozen observer files unchanged',
        environment_scope=environment['scope'],validated_environment_blocks=environment['validated_blocks'],
        source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in sources+['audit_awake_environment.py']}),indent=2)+'\n')
