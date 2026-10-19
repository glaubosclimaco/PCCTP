"""Passive clock-duration observations; no causal attribution of prior stalls."""
import time,json,datetime
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    target=HERE/'clock_monitor_observations.json';assert not target.exists(),'Preserve previous observations; do not duplicate'
    rows=[]
    for i in range(6):
        native=time.clock_gettime(time.CLOCK_MONOTONIC);python=time.monotonic();perf=time.perf_counter()
        row=dict(index=i,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),native_clock=native,python_monotonic=python,perf_counter=perf,parent_cpu=time.process_time())
        if rows:
            previous=rows[-1];row.update(native_interval=native-previous['native_clock'],python_interval=python-previous['python_monotonic'],perf_interval=perf-previous['perf_counter'])
            row['duration_difference_seconds']=row['native_interval']-row['python_interval']
        rows.append(row);target.write_text(json.dumps(dict(scope='passive operational monitor; no algorithm benchmark',complete=i==5,rows=rows),indent=2)+'\n')
        if i<5:time.sleep(30)
    print(json.dumps(dict(samples=len(rows),maximum_absolute_duration_difference=max(abs(r.get('duration_difference_seconds',0)) for r in rows))))
