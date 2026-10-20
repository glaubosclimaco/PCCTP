"""Compare recorded public metadata checksums with approved TRAIN files only."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
META = ROOT/'research_notes/benchmark_metadata/mendeley_train_identity.json'

if __name__ == '__main__':
    metadata = json.loads(META.read_text())
    profile = json.loads((ROOT/'experiments/lkh/hard_train_instance_profile.json').read_text())
    allowed = {row['instance']: row for row in profile['rows']}
    roles_path = ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json'
    roles = json.loads(roles_path.read_text())['roles']
    assert len(allowed) == 9
    assert all(roles[row['geometry']] == 'train' for row in allowed.values())
    seen = set()
    for row in metadata['rows']:
        name = row['instance']
        assert name in allowed and name not in seen
        seen.add(name)
        assert row['source_url'] == ('https://data.mendeley.com/datasets/8yjbfgcfvn/2/files/' + row['file_id'])
        sha = hashlib.sha256((ROOT/'data'/name).read_bytes()).hexdigest()
        assert sha == row['public_sha256'] == row['local_sha256'] == allowed[name]['input_sha256']
        assert row['identical'] is True
    assert metadata['unverified_remaining_train_instances'] == 9-len(seen)
    report = dict(status='COMPLETE_TRAIN' if len(seen) == 9 else 'PARTIAL_TRAIN',
                  expected_train_instances=9, verified_train_instances=len(seen),
                  unverified_train_instances=sorted(set(allowed)-seen),
                  scope='Public SHA-256 metadata matches only; no reserved data, downloads or solver calls. '
                        'Not a reproduction of the published 107-case evaluation.',
                  metadata_sha256=hashlib.sha256(META.read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    META.with_name('mendeley_train_identity_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
