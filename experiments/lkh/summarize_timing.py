"""Describe corrected timing diagnostics without assigning unmeasured causes."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    folder=HERE/'timing_diagnostic_v2';assert (folder/'completion.json').exists()
    rows=[]
    for p in sorted(folder.glob('g_*.json')):
        if p.name.endswith('.native.json'):continue
        r=json.loads(p.read_text());n=r['native'];cpu=r['parent_cpu_seconds']+r['child_user_cpu_seconds']+r['child_system_cpu_seconds']
        rows.append(dict(instance=r['instance'],wall_seconds=r['wall_seconds'],cpu_seconds=cpu,cpu_wall_ratio=cpu/r['wall_seconds'],startup_debit_seconds=n['startup_debit_seconds'],native_incumbent=n['feasible_before_deadline'],native_cost=n['cost'],hard_timeout=n['hard_timeout'],calibration_wall_cpu_before=r['calibration_before']['wall_seconds']/r['calibration_before']['cpu_seconds'],calibration_wall_cpu_after=r['calibration_after']['wall_seconds']/r['calibration_after']['cpu_seconds']))
    report=dict(scope='three corrected operational observations, not algorithm benchmark or explanation of prior stalls',rows=rows,interpretation='CPU and elapsed clock describe these runs only. Prior batches lacked CPU observations; no retrospective attribution to load, sleep, or app behavior is justified.')
    (HERE/'timing_diagnostic_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
