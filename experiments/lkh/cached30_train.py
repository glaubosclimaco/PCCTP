"""Prospective conditional cache collection; blocked until priority studies finish.

This module is not queued by any active worker. Running it acquires the same
kernel lock as paired70, verifies prerequisites and certificates before output
creation, and preserves interrupted pairs for review.
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from experiments.lkh.cached30_starts import load_starts
from experiments.lkh.cached30_phase import source_hashes
from experiments.lkh.cached30_pairs import execute_pairs
from experiments.lkh.cached30_batch_audit import audit_batch
from experiments.lkh.paired70_train import study_lock

OUT = HERE/'cached30_train'
CERTIFICATES = ('verify_cached_polish.json', 'verify_cached30_phase.json',
                'verify_cached30_pairs.json', 'verify_cached30_batch.json', 'verify_cached30_owner.json')


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_certificates():
    for name in CERTIFICATES:
        certificate = json.loads((HERE/name).read_text())
        assert certificate['status'] == 'PASS'
        for filename, hashed in certificate['source_sha256'].items():
            assert digest(ROOT/filename) == hashed, 'Stale certificate: '+name+' / '+filename
    lock = json.loads((HERE/'verify_paired70_lock.json').read_text())
    assert lock['status'] == 'PASS' and lock['source_sha256'] == digest(HERE/'paired70_train.py')


def verify_prerequisites():
    # Do not run even read-only full auditors until both priority batches ended.
    for folder in ('ruin_total300_awake_train', 'published70_calibration'):
        assert (HERE/folder/'completion.json').exists(), 'Wait for completed '+folder
    subprocess.run([sys.executable, str(HERE/'observe_awake300.py')], cwd=ROOT, check=True)
    subprocess.run([sys.executable, str(HERE/'audit_published70_calibration.py')], cwd=ROOT, check=True)
    paths = ['ruin_total300_awake_train_audit.json', 'ruin_total300_awake_environment_audit.json',
             'published70_calibration_audit.json']
    awake, environment, calibration = [json.loads((HERE/name).read_text()) for name in paths]
    assert awake['status'] == 'PASS' and awake['batch_complete'] and awake['runs'] == 81
    assert environment['batch_complete'] and environment['eligible_complete_environment']
    assert environment['validated_blocks'] == environment['expected_blocks'] == 27
    assert calibration['status'] == 'PASS' and calibration['batch_complete']
    assert calibration['runs'] == calibration['completed70'] == 9 and calibration['all_references_completed']
    paths += ['ruin_total300_awake_train/protocol.json', 'ruin_total300_awake_train/completion.json',
              'published70_calibration/protocol.json', 'published70_calibration/completion.json']
    return {name: digest(HERE/name) for name in paths}


def verify_ordering_binary(binary):
    # Calibration runs MDM only: its protocol intentionally has no LKH entry.
    awake = json.loads((HERE/'ruin_total300_awake_train/protocol.json').read_text())
    assert awake['binary_hashes'].get(str(binary)) == digest(binary), 'Ordering binary changed since awake300'


def collect():
    prerequisites = verify_prerequisites()
    verify_certificates()
    manifest = load_starts()
    assert json.loads((HERE/'cached30_starts.json').read_text()) == manifest
    starts = manifest['starts']
    phase_sources = source_hashes()
    files = set(phase_sources) | {'experiments/lkh/'+name for name in (
        'cached30_train.py', 'cached30_pairs.py', 'cached30_starts.py', 'cached30_starts.json',
        'cached30_phase_audit.py', 'cached30_batch_audit.py', 'raw_route_audit.py',
        'awake300_train.py', 'audit_awake_environment.py', 'paired70_train.py')}
    files |= {'experiments/lkh/'+name for name in CERTIFICATES}
    files.add('experiments/lkh/verify_paired70_lock.json')
    files.add('research_notes/cached30_conditional_design.md')
    binary = Path('/private/tmp/pcctp-lkh3-research/LKH')
    verify_ordering_binary(binary)
    protocol = dict(scope=manifest['scope'], starts=starts, budget_seconds=30.,
                    scheduled_pairs=27, scheduled_phases=54, arms=['uncached', 'cached'],
                    repair_seed_rule='source_seed+49', prerequisite_sha256=prerequisites,
                    source_sha256={name: digest(ROOT/name) for name in sorted(files)},
                    phase_source_sha256=phase_sources, binary_sha256={str(binary): digest(binary)})
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT/'protocol.json'
    if path.exists(): assert json.loads(path.read_text()) == protocol, 'Frozen protocol changed'
    else: path.write_text(json.dumps(protocol, indent=2)+'\n')
    protocol_hash = digest(path)
    # Independently check all committed pairs before allowing resume.
    audit_batch(OUT, starts, protocol_hash, phase_sources, require_complete=False)
    execute_pairs(starts, OUT, protocol_hash, phase_sources)
    report = audit_batch(OUT, starts, protocol_hash, phase_sources)
    audit_path = HERE/'cached30_train_audit.json'
    audit_path.write_text(json.dumps(report, indent=2)+'\n')
    completion = dict(scheduled_pairs=27, scheduled_phases=54, audited_pairs=report['audited_pairs'],
                      audited_phases=report['audited_phases'], protocol_sha256=protocol_hash,
                      audit_sha256=digest(audit_path))
    temp = OUT/'completion.tmp'; temp.write_text(json.dumps(completion, indent=2)+'\n')
    temp.replace(OUT/'completion.json')
    print(json.dumps({k: report[k] for k in ('status', 'audited_pairs', 'mean_gain_percent',
                                            'wins', 'losses', 'ties', 'gate_passed')}), flush=True)


def main():
    with study_lock(HERE/'paired70_train.queue.lock'):
        collect()


if __name__ == '__main__': main()
