"""Single durable queue: audit completed awake300, then run seed80 calibration."""
import fcntl
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def readiness():
    out = HERE/'ruin_total300_awake_train'
    failed = [p.name for p in out.glob('environment.g_*.json')
              if json.loads(p.read_text())['status'] != 'PASS']
    return dict(awake_completed=(out/'completion.json').exists(), failed_blocks=failed,
                calibration_completed=(HERE/'published70_calibration/completion.json').exists())


def main():
    if sys.argv[1:] == ['--check']:
        print(json.dumps(readiness())); return
    assert len(sys.argv) == 1
    with (HERE/'published70_awake_calibration.queue.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        lock.seek(0); lock.truncate(); lock.write(str(os.getpid())); lock.flush()
        print(json.dumps(dict(status='WAITING_FOR_AWAKE81', pid=os.getpid())), flush=True)
        while True:
            state = readiness()
            assert not state['failed_blocks'], 'Review compromised awake block before restarting queue'
            if state['awake_completed']: break
            time.sleep(20)
        subprocess.run([sys.executable, str(HERE/'observe_awake300.py')], cwd=ROOT, check=True)
        env = json.loads((HERE/'ruin_total300_awake_environment_audit.json').read_text())
        audit = json.loads((HERE/'ruin_total300_awake_train_audit.json').read_text())
        assert env['eligible_complete_environment'] and env['validated_blocks'] == 27
        assert audit['status'] == 'PASS' and audit['batch_complete'] and audit['runs'] == 81
        if not state['calibration_completed']:
            subprocess.run(['/usr/bin/caffeinate', '-i', sys.executable,
                            str(HERE/'published70_calibration.py')], cwd=ROOT, check=True)
        subprocess.run([sys.executable, str(HERE/'audit_published70_calibration.py')], cwd=ROOT, check=True)
        print(json.dumps(dict(status='CALIBRATION_FINISHED_AND_AUDITED')), flush=True)


if __name__ == '__main__': main()
