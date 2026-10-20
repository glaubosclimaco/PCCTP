"""Real-process lock contention and exception release; no optimizers or inputs."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from paired70_train import study_lock

if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    code = ('import sys; from pathlib import Path; '
            'from paired70_train import study_lock; '
            '\nwith study_lock(Path(sys.argv[1])): pass')
    with tempfile.TemporaryDirectory(prefix='pcctp-paired-lock-', dir='/private/tmp') as tmp:
        path = Path(tmp) / 'queue.lock'
        with study_lock(path):
            original = path.read_bytes()
            denied = subprocess.run([sys.executable, '-c', code, str(path)], cwd=here,
                                    capture_output=True, text=True, timeout=10)
            assert denied.returncode != 0 and 'already running' in denied.stderr
            assert path.read_bytes() == original, 'Contender changed owner metadata'
        allowed = subprocess.run([sys.executable, '-c', code, str(path)], cwd=here,
                                 capture_output=True, text=True, timeout=10)
        assert allowed.returncode == 0, allowed.stderr
        try:
            with study_lock(path):
                raise ValueError('synthetic collection failure')
        except ValueError:
            pass
        with study_lock(path):
            assert json.loads(path.read_text())['pid'] > 0
        assert path.exists(), 'Lock inode must remain for safe reuse'
    report = dict(status='PASS', scope='Synthetic process locking only; no optimizer or benchmark input',
                  checks=['second process rejected before entering critical section',
                          'rejected process leaves owner metadata unchanged',
                          'new process acquires lock after normal release',
                          'exception releases lock', 'lock inode retained'],
                  source_sha256=hashlib.sha256((here/'paired70_train.py').read_bytes()).hexdigest())
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
