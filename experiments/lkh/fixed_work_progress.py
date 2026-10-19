"""Describe partial fixed-work results without averaging censored and completed runs."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    audit=json.loads((HERE/'fixed_work_audit.json').read_text());assert audit['status']=='PASS'
    rows=audit['rows'];completed=[r for r in rows if r['status']=='COMPLETED70'];censored=[r for r in rows if r['status']=='CENSORED'];missing=9-len(rows)
    lines=['# Fixed-work diagnostic: progress, not final results','',f'{len(completed)} completed 70 iterations; {len(censored)} censored; {missing} pending. No aggregate performance conclusion is drawn from this partial set.','','| Geometry | Status | Last saved iteration | MDM cost | Additional LKH cost | CPU (min) | Wall (min) | Supplied MDM reference |','|---|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['status']=='NO_VALID_EXPORT':continue
        lines.append(f"| {r['instance'].split('_VT')[0]} | {r['status']} | {r['last_iteration']} | {r['mdm_cost']:g} | {r['lkh_cost']:g} | {r['cpu_seconds']/60:.2f} | {r['native_wall_seconds']/60:.2f} | {r['published_best']:g} |")
    lines+=['','The last saved iteration of a censored run is a checkpoint, not an assertion about the exact iteration at interruption. The requested external timeout was 1800 seconds; observed return times are reported separately and can exceed it. LKH received up to two additional seconds, so its cost change is not an equal-total-time comparison. Published reference costs are not assumed to be certified optima.']
    (ROOT/'research_notes/FIXED_WORK_PROGRESS.md').write_text('\n'.join(lines)+'\n')
    print(f'Progress recorded: {len(completed)} complete, {len(censored)} censored, {missing} pending')
