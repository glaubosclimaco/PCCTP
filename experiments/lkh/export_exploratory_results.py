"""Export only complete audited studies; no active or reserved results read."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

if __name__=='__main__':
 names=['ruin_total60_absolute_summary.json','ruin_total60_absolute_audit.json',
        'fixed_support_ablation_audit.json','ruin_validation60_summary.json',
        'ruin_validation60_audit.json','joint_insertion_train_summary.json',
        'joint_insertion_train_audit.json']
 data={n:json.loads((HERE/n).read_text()) for n in names}
 for n,d in data.items():
  if n.endswith('_audit.json'):assert d['status']=='PASS' and d['batch_complete']
  else:assert d['status']=='COMPLETE'
 train=data[names[0]];fixed=data[names[2]];val=data[names[3]];joint=data[names[5]]
 fg=[r['hybrid_gain_percent'] for r in fixed['geometry']]
 jg=[r['gain_vs_greedy'] for r in joint['rows']]
 rows=[('Charged TRAIN, 60 s',54,train['mean_gain_percent'],train['wins'],train['losses'],9),
       ('Conditional fixed-support replay',27,sum(fg)/len(fg),sum(g>1e-9 for g in fg),sum(g< -1e-9 for g in fg),9),
       ('Limited validation, 60 s',70,val['mean_gain_percent'],val['wins'],val['losses'],7),
       ('Conditional joint insertion vs. greedy',324,joint['mean_gain_vs_greedy_percent'],joint['wins_vs_greedy'],sum(g< -1e-9 for g in jg),9)]
 lines=[r'\begin{table}[htbp]',r'\centering',r'\small',
        r'\begin{tabular}{lrrr}',r'\toprule',
        r'Study & Runs & Reduction (\%) & W/L/T \\',r'\midrule']
 for label,runs,effect,w,l,n in rows:lines.append(f'{label} & {runs} & {effect:.3f} & {w}/{l}/{n-w-l} '+r'\\')
 lines += [r'\bottomrule',r'\end{tabular}',
 r'\caption{Completed exploratory comparisons. Positive reduction favors the hybrid in the first three rows and nearest-gap joint insertion in the last. W/L/T counts geometry means. Conditional rows exclude generation of their common starting states and are not directly comparable with charged end-to-end rows.}',
 r'\label{tab:exploratory-complete}',r'\end{table}']
 target=ROOT/'submission/tables/exploratory_complete.tex';target.parent.mkdir(exist_ok=True)
 target.write_text('\n'.join(lines)+'\n')
 provenance=dict(scope='complete exploratory studies only; no active300s/test results',
  source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
  output_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
  rows=[dict(study=x[0],runs=x[1],effect=x[2],wins=x[3],losses=x[4],ties=x[5]-x[3]-x[4]) for x in rows])
 target.with_suffix('.provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 print(json.dumps(dict(status='PASS',rows=len(rows),output=str(target))))
