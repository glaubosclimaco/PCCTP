"""One common-scale figure for both arms; completed audited TRAIN data only."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def collect_plot_data():
    closing_path = HERE/'awake300_closing_analysis.json'
    assert closing_path.exists(), 'Wait for completed closing audit before rendering'
    closing = json.loads(closing_path.read_text())
    assert closing['status'] == 'COMPLETE_DESCRIPTIVE_TRAIN' and closing['runs'] == 81
    audit_path = HERE/'ruin_total300_awake_train_audit.json'
    assert digest(audit_path) == closing['upstream_sha256'][str(audit_path.relative_to(ROOT))]
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'PASS' and audit['batch_complete'] and len(audit['rows']) == 81
    records = {(r['instance'], r['seed'], r['mode']): r for r in audit['rows']}
    assert len(records) == 81 and len(closing['geometries']) == 9
    data = []
    for geometry in closing['geometries']:
        name = geometry['instance']; means, points = {}, {}
        costs = {mode: [Fraction(str(records[name, seed, mode]['cost'])) for seed in (70,71,72)]
                 for mode in ('mdm','hybrid10','hybrid60')}
        for mode in ('hybrid10','hybrid60'):
            gain = 100*(sum(costs['mdm'])-sum(costs[mode]))/sum(costs['mdm'])
            assert gain == Fraction(geometry['exact_gains'][mode])
            means[mode] = float(gain)
            points[mode] = [float(100*(a-b)/a) for a,b in zip(costs['mdm'],costs[mode])]
        data.append(dict(instance=name, means=means, paired_seed_effects=points))
    return data, closing_path, audit_path


def main():
    data, closing_path, audit_path = collect_plot_data()
    # Avoid creating any rendering output or cache before complete-data gates.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    colors = {'hybrid10': '#D55E00', 'hybrid60': '#0072B2'}
    offsets = {'hybrid10': -.18, 'hybrid60': .18}
    allvalues = [x for row in data for values in row['paired_seed_effects'].values() for x in values]
    allvalues += [x for row in data for x in row['means'].values()]
    lo, hi = min(0., min(allvalues)), max(0., max(allvalues))
    margin = max(1., .08*(hi-lo))
    ax.set_xlim(lo-margin, hi+margin)
    for i, row in enumerate(data):
        for mode in ('hybrid10','hybrid60'):
            y, gain = i+offsets[mode], row['means'][mode]
            ax.barh(y, gain, height=.26, color=colors[mode], alpha=.8, zorder=2)
            for j, marker in enumerate(('o','s','D')):
                ax.scatter(row['paired_seed_effects'][mode][j], y+(j-1)*.055, marker=marker,
                           s=15, c='#222222', alpha=.55, edgecolors='white', linewidths=.35, zorder=3)
            ax.annotate(f'{gain:+.2f}', (gain,y), xytext=(4 if gain>=0 else -4,0),
                        textcoords='offset points', ha='left' if gain>=0 else 'right',
                        va='center', fontsize=8, color=colors[mode], zorder=4)
    ax.set_yticks(range(9), [r['instance'].split('_VT')[0].removeprefix('g_') for r in data])
    ax.invert_yaxis(); ax.axvline(0, color='#666666', linewidth=.8)
    ax.grid(axis='x', color='#e1e1e1', linewidth=.6); ax.set_axisbelow(True)
    for spine in ('top','right','left'): ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('Cost reduction versus paired MDM300 (%)')
    ax.set_title('TRAIN only: both prespecified hybrids, 300 s total', fontsize=12, pad=14)
    handles = [Line2D([0],[0], color=colors[m], linewidth=6, label=label)
               for m,label in (('hybrid10','MDM10 + repair'),('hybrid60','MDM60 + repair'))]
    handles += [Line2D([0],[0], color='#333333', marker=marker, linestyle='None', markersize=4,
                      label=f'Seed {70+j}') for j,marker in enumerate(('o','s','D'))]
    fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=8,
               bbox_to_anchor=(.54,.025))
    fig.subplots_adjust(left=.13,right=.97,top=.9,bottom=.18)
    out = ROOT/'submission/figures'; out.mkdir(parents=True,exist_ok=True)
    for extension in ('png','svg'): fig.savefig(out/f'awake300_effects.{extension}',dpi=240)
    caption = ('Selected TRAIN geometries, three paired seeds, 300 seconds including initialization. '
               'Bars show reductions of seed-mean costs; markers show each seed-pair reduction, '
               'not confidence intervals. One shared axis includes all losses and both arms. '
               'Positive values favor the hybrid. This is an exploratory short-horizon screen, '
               'not independent test evidence or superiority over completed MDM70.')
    provenance = dict(scope=caption, data=data, source_sha256={str(p.relative_to(ROOT)):digest(p)
                      for p in (Path(__file__),closing_path,audit_path)},
                      output_sha256={f'awake300_effects.{ext}':digest(out/f'awake300_effects.{ext}')
                                    for ext in ('png','svg')})
    (out/'awake300_effects.provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    plt.close(fig)
    print('Saved complete TRAIN figure:18mean bars/54paired-seed markers; visual inspection still required')


if __name__ == '__main__': main()
