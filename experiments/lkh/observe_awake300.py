"""Reuse unchanged original observers with explicit fresh-batch output paths."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
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
    (HERE/'awake300_observer_provenance.json').write_text(json.dumps(dict(
        transformation='literal batch/output-name substitution only; original frozen observer files unchanged',
        source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in sources}),indent=2)+'\n')
