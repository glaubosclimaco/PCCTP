"""Observed checkpoints only: no interpolated anytime trajectory."""
import os,json,hashlib
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/pcctp-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
if __name__=='__main__':
    name='g_a280_VT56_T168_W56_50';folder=ROOT/'experiments/lkh/fixed_work';p=folder/f'{name}.json';r=json.loads(p.read_text());assert r['status']=='COMPLETED70'
    raw=folder/f'{name}.native.raw';payload=json.loads(raw.read_text());points=payload['snapshots'];reference=r['published_mdm_best']
    fig,ax=plt.subplots(figsize=(7.3,4.2));x=[s['seconds']/60 for s in points];y=[100*(s['cost']-reference)/reference for s in points]
    ax.scatter(x,y,color='#0072B2',s=48,zorder=3)
    for i,s in enumerate(points):ax.annotate(f"{s['iteration']} iter.",(x[i],y[i]),xytext=(0,9 if i!=4 else 18),textcoords='offset points',ha='center',fontsize=9)
    ax.axhline(0,color='#666666',linewidth=1,linestyle='--')
    ax.set(title='a280: the published MDM cost is reached by checkpoint 20',xlabel='Elapsed wall time (minutes)',ylabel='Gap to supplied published cost (%)')
    ax.set_ylim(-.4,max(y)+1.3);ax.set_xlim(0,max(x)+2);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.grid(axis='y',alpha=.2)
    fig.text(.12,.02,'Observed checkpoints only; no interpolation. One training instance, seed 42.\nReference cost 2001 is not an optimality certificate. LKH after iteration 70: no improvement.',fontsize=8,color='#444444')
    fig.tight_layout(rect=(0,.1,1,1));out=ROOT/'submission/figures';out.mkdir(exist_ok=True,parents=True)
    fig.savefig(out/'a280_fixed_work_checkpoints.png',dpi=180);fig.savefig(out/'a280_fixed_work_checkpoints.svg');plt.close(fig)
    (out/'a280_fixed_work_checkpoints.provenance.json').write_text(json.dumps(dict(scope='diagnostic figure, not final paper evidence',input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),reference=reference,points=points),indent=2)+'\n')
