"""Read-only SHA256 check of the declared completed-study artifact scope."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    hashed = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''): hashed.update(chunk)
    return hashed.hexdigest()


if __name__ == '__main__':
    manifest = json.loads((ROOT/'submission/reproducibility_manifest.json').read_text())
    issues = []
    for filename, expected in manifest['artifact_sha256'].items():
        path = (ROOT/filename).resolve()
        path.relative_to(ROOT)
        if not path.is_file(): issues.append(dict(path=filename,status='MISSING'))
        elif digest(path) != expected: issues.append(dict(path=filename,status='CHANGED'))
    runtime = []
    for filename, expected in manifest['external_binary_sha256'].items():
        path = Path(filename)
        status = 'MISSING' if not path.is_file() else ('PASS' if digest(path)==expected else 'CHANGED')
        runtime.append(dict(path=filename,status=status))
    print(json.dumps(dict(artifact_status='PASS' if not issues else 'FAIL',
        artifacts_checked=len(manifest['artifact_sha256']),artifact_issues=issues,
        external_runtime=runtime,scope=manifest['scope']),indent=2))
    if issues: raise SystemExit(1)
