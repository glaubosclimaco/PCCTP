"""Audit only the separately restarted study; preserve all original reports."""
import hashlib,json
from pathlib import Path
from audit_paired70_train import audit_batch
from summarize_paired70_train import summarize
from paired70_dispersion import report as dispersion
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    folder=HERE/'paired70_rebuilt_train';audit=audit_batch(folder)
    if audit['status']=='NO_RESULTS':print('NO_RESULTS');raise SystemExit(0)
    protocol=json.loads((folder/'protocol.json').read_text())
    assert protocol['superseded_incomplete_study']=='paired70_train'
    assert protocol['binary_hashes']['/private/tmp/pcctp-lkh3-research/LKH']=='1908f80deeafc8952acb21d36ba47ed618a71f2c4527c9654338cd24316be949'
    audit['scope']+='; separate rebuilt-LKH restart, old five pairs excluded'
    path=HERE/'paired70_rebuilt_audit.json';path.write_text(json.dumps(audit,indent=2)+'\n')
    summary=summarize(audit,protocol);summary['audit_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    (HERE/'paired70_rebuilt_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (HERE/'paired70_rebuilt_dispersion.json').write_text(json.dumps(dispersion(audit,protocol),indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='rows'}))
