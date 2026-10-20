"""Training-only paired effects; never reads validation or test outcomes."""
import json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
if __name__=='__main__':
 a=HERE/'ruin_total60_absolute_summary.json';b=HERE/'fixed_support_ablation_audit.json'
 total=json.loads(a.read_text());ablation=json.loads(b.read_text());assert total['status']=='COMPLETE' and ablation['batch_complete'] and ablation['status']=='PASS'
 names=[r['instance'] for r in total['rows']];assert len(names)==9 and set(names)=={r['instance'] for r in ablation['geometry']}
 second={r['instance']:r for r in ablation['geometry']};pairs={(r['instance'],r['seed']):r for r in ablation['rows']}
 means=[[],[]];points=[[],[]];data=[]
 for row in total['rows']:
  name=row['instance'];one=[100*(p['mdm_cost']-p['hybrid_cost'])/p['mdm_cost'] for p in row['pairs']];two=[]
  for seed in (42,43,44):
   p=pairs[(name,seed)];two.append(100*(p['fixed_cost']-p['hybrid_cost'])/p['fixed_cost'])
  means[0].append(row['gain_percent']);means[1].append(second[name]['hybrid_gain_percent']);points[0].append(one);points[1].append(two)
  data.append(dict(instance=name,total_gain=row['gain_percent'],conditional_gain=second[name]['hybrid_gain_percent'],seed_total_gains=one,seed_conditional_gains=two))
 plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,axes=plt.subplots(1,2,figsize=(9.0,4.8),sharey=True)
 allvalues=[v for panel in points for row in panel for v in row]+[v for panel in means for v in panel];low=min(-.5,min(allvalues)-.3);high=max(allvalues)+.5
 for j,ax in enumerate(axes):
  y=np.arange(9);ax.barh(y,means[j],height=.55,color=['#0072B2' if v>=0 else '#D55E00' for v in means[j]],alpha=.8,zorder=2)
  for seed_index,marker in enumerate(('o','s','D')):
   ax.scatter([p[seed_index] for p in points[j]],y+(seed_index-1)*.12,s=19,marker=marker,c='#202020',edgecolors='white',linewidths=.4,zorder=3,label=f'Seed {42+seed_index}')
  ax.axvline(0,color='#777777',linewidth=.8);ax.grid(axis='x',color='#e2e2e2',linewidth=.6,zorder=0);ax.set_axisbelow(True);ax.set_xlim(low,high)
  for spine in ('top','right','left'):ax.spines[spine].set_visible(False)
  ax.tick_params(axis='y',length=0);ax.set_xlabel('Hybrid cost reduction (%)');ax.set_title(('A  Versus MDM, 60 s total','B  Versus fixed-support LKH replay')[j],fontsize=10,pad=12)
 axes[0].set_yticks(np.arange(9),[name.split('_VT')[0].removeprefix('g_') for name in names]);axes[0].invert_yaxis()
 handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=3,loc='lower center',bbox_to_anchor=(.52,.01),frameon=False,fontsize=9)
 fig.suptitle('Training only: nine geometries, three paired seeds',fontsize=12,y=.98);fig.subplots_adjust(left=.11,right=.98,top=.86,bottom=.19,wspace=.14)
 out=ROOT/'submission/figures';out.mkdir(exist_ok=True)
 for extension in ('png','svg','pdf'):fig.savefig(out/f'hybrid_training_effects.{extension}',dpi=220)
 caption='Training screen on nine geometries and three seeds. Bars: relative reduction of seed-mean costs; markers: individual seed-pair reductions, not confidence intervals. A: fresh end-to-end comparison with MDM under60seconds including initialization. B: conditional fixed-support LKH replay from the exact saved pre-repair starts and remaining allowances. These panels have different controls; neither is independent test evidence.'
 provenance=dict(scope='training only; no validation/test files read',caption=caption,source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (a,b,Path(__file__))},data=data)
 (out/'hybrid_training_effects.provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');print('Saved training-only effects figure,9geometries x3seeds')
