"""Synthetic exposure-bound fixtures, independent of all benchmark data."""
import hashlib,json,tempfile
from pathlib import Path
from profile_mining_exposure import profile

if __name__=='__main__':
 with tempfile.TemporaryDirectory(prefix='pcctp-mining-bounds-',dir='/private/tmp') as tmp:
  folder=Path(tmp)
  cases=[('completed_overrun',2,12.,True,1,8.,1,2),
         ('checkpoint_only_overrun',None,None,False,2,12.,0,None),
         ('completed_no_calls',0,12.,True,0,9.,0,0)]
  for index,(label,calls,seconds,complete,cc,cs,lower,upper) in enumerate(cases):
   stem=f'g_synthetic{index}.seed1.mdm';raw=folder/(stem+'.native.raw');checkpoint=Path(str(raw)+'.checkpoint')
   if calls is not None:raw.write_text(json.dumps(dict(mining_calls=calls,seconds=seconds,completed=complete,iterations=10)))
   else:raw.write_text('{truncated')
   checkpoint.write_text(json.dumps(dict(mining_calls=cc,seconds=cs,completed=False,iterations=5)))
   selected=raw if calls is not None else checkpoint
   row=dict(instance=label,seed=1,mode='mdm',native=dict(budget_seconds=10,payload_source=selected.name))
   (folder/(stem+'.json')).write_text(json.dumps(row))
   (folder/(stem+'.native.log')).write_text('Minerador executado.\n'*(calls if calls is not None else cc))
  report=profile(folder)
  for row,case in zip(report['rows'],cases):
   assert row['instance']==case[0]
   assert row['timely_mining_calls_lower_bound']==case[6]
   assert row['timely_mining_calls_upper_bound']==case[7]
  assert len(report['rows'][1]['payload_read_errors'])==1
  assert report['modes']['mdm']['runs_with_proven_timely_call']==1
  assert report['modes']['mdm']['runs_with_unknown_upper_bound']==1
 source=Path(__file__).with_name('profile_mining_exposure.py')
 output=dict(status='PASS',synthetic_cases=3,scope='finite bound/partial-export checks only',
             profiler_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
 Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output))
