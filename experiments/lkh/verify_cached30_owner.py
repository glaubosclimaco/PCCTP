"""Owner lock-route and prerequisite checks; stubs only, no collection."""
import hashlib
import json
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.lkh import cached30_train as module


def main():
    original = module.HERE, module.ROOT, module.OUT, module.study_lock, module.collect, module.subprocess.run
    checks, audit_calls = [], []
    state = dict(held=False, path=None, visited=0)
    @contextmanager
    def lock(path):
        state.update(held=True, path=path)
        try: yield
        finally: state['held'] = False
    def collected():
        assert state['held']
        state['visited'] += 1
    def rejected(function):
        try: function()
        except AssertionError: pass
        else: raise AssertionError('Expected prerequisite/certificate rejection')
    def write(path, value): path.write_text(json.dumps(value)+'\n')
    try:
        with tempfile.TemporaryDirectory(prefix='pcctp-cache-owner-') as folder:
            root = Path(folder); here = root/'experiments/lkh'; here.mkdir(parents=True)
            module.ROOT, module.HERE, module.OUT = root, here, here/'cached30_train'
            module.study_lock, module.collect = lock, collected
            module.main()
            assert state['path'] == here/'paired70_train.queue.lock' and state['visited'] == 1 and not state['held']
            checks.append('collection entered only while shared paired70 lock context held')
            def interrupted():
                assert state['held']
                raise RuntimeError('synthetic owner interruption')
            module.collect = interrupted
            try: module.main()
            except RuntimeError: pass
            else: raise AssertionError('Owner exception swallowed')
            assert not state['held']
            checks.append('owner exception exits shared lock context')
            module.collect = original[4]
            module.subprocess.run = lambda *a, **k: audit_calls.append(a[0])
            rejected(module.collect)
            assert not audit_calls and not module.OUT.exists()
            checks.append('missing awake completion blocks before auditors, outputs or optimizer')
            awake = here/'ruin_total300_awake_train'; awake.mkdir()
            write(awake/'completion.json', {}); write(awake/'protocol.json', {})
            rejected(module.collect)
            assert not audit_calls and not module.OUT.exists()
            checks.append('missing calibration completion blocks before auditors, outputs or optimizer')
            calibration = here/'published70_calibration'; calibration.mkdir()
            write(calibration/'completion.json', {}); write(calibration/'protocol.json', {})
            aa = here/'ruin_total300_awake_train_audit.json'
            ea = here/'ruin_total300_awake_environment_audit.json'
            ca = here/'published70_calibration_audit.json'
            a = dict(status='PASS', batch_complete=True, runs=81)
            e = dict(batch_complete=True, eligible_complete_environment=True, validated_blocks=27, expected_blocks=27)
            c = dict(status='PASS', batch_complete=True, runs=9, completed70=9, all_references_completed=True)
            for path, payload in ((aa, a), (ea, e), (ca, c)): write(path, payload)
            hashes = module.verify_prerequisites()
            assert len(audit_calls) == 2 and len(hashes) == 7
            assert audit_calls[0][-1].endswith('observe_awake300.py')
            assert audit_calls[1][-1].endswith('audit_published70_calibration.py')
            checks.append('both completions cause fresh read-only audits and bind seven provenance files')
            write(aa, dict(a, runs=80)); rejected(module.verify_prerequisites); write(aa, a)
            checks.append('80 of81 awake results cannot satisfy launch prerequisite')
            write(ea, dict(e, eligible_complete_environment=False)); rejected(module.verify_prerequisites); write(ea, e)
            checks.append('compromised awake environment blocks launch')
            write(ca, dict(c, completed70=8, all_references_completed=False)); rejected(module.verify_prerequisites); write(ca, c)
            checks.append('censored or incomplete full reference blocks launch')
            source = root/'synthetic.py'; source.write_text('synthetic source\n')
            for name in module.CERTIFICATES:
                write(here/name, dict(status='PASS', source_sha256={'synthetic.py': module.digest(source)}))
            pair_source = here/'paired70_train.py'; pair_source.write_text('synthetic lock source\n')
            write(here/'verify_paired70_lock.json', dict(status='PASS', source_sha256=module.digest(pair_source)))
            module.verify_certificates()
            checks.append('current source-bound certificates accepted')
            source.write_text('changed\n'); rejected(module.verify_certificates)
            checks.append('changed adapter/owner certificate source blocks launch')
            source.write_text('synthetic source\n'); pair_source.write_text('changed lock\n')
            rejected(module.verify_certificates)
            checks.append('changed shared kernel-lock implementation blocks launch')
    finally:
        module.HERE, module.ROOT, module.OUT, module.study_lock, module.collect, module.subprocess.run = original
    files = ['experiments/lkh/'+name for name in ('cached30_train.py', 'verify_cached30_owner.py',
             'cached30_pairs.py', 'cached30_starts.py', 'cached30_batch_audit.py',
             'cached30_phase.py', 'cached30_phase_audit.py', 'paired70_train.py')]
    report = dict(status='PASS', scope='Synthetic owner lock routing, launch prerequisites and source '
                  'certificates; no benchmark inputs, real audit subprocesses or optimizer calls. '
                  'Kernel exclusion itself is covered by verify_paired70_lock.json.', checks=checks,
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', checks=len(checks), native_optimizer_calls=0)))


if __name__ == '__main__': main()
