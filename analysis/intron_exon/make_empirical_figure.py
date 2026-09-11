"""Plot the previously completed direct paired-bootstrap profile offsets."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent

def main():
    df = pd.read_csv(ROOT / 'empirical' / 'stage_lag_bootstrap.csv')
    order = [('KR', 1), ('MLPT', 1), ('MLPT', 2), ('SVB', 1), ('SVB', 2)]
    labels = ['Kr', 'mlpt, pulse 1', 'mlpt, pulse 2', 'svb, pulse 1', 'svb, pulse 2']
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.linewidth': .7})
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), sharey=True)
    for ax, metric, title, color in zip(axes, ['rising', 'falling'],
            ['Rising half-height', 'Falling half-height'], ['#256a9a', '#a74b39']):
        ax.axvline(0, color='.45', lw=1, ls='--')
        for i, (gene, pulse) in enumerate(order):
            row = df[(df.gene == gene) & (df.pulse == pulse) & (df.metric == metric)].iloc[0]
            ax.errorbar(row.estimate, i, xerr=[[row.estimate-row.lo95], [row.hi95-row.estimate]],
                        fmt='o', color=color, capsize=3, lw=1.4, ms=5)
        ax.set_title(title, loc='left', weight='bold', fontsize=11)
        ax.set_xlim(-1.2, 1.7)
        ax.set_xlabel('Exon minus intron crossing time\n(eve-en substage units)')
        ax.grid(axis='x', alpha=.12)
    axes[0].set_yticks(np.arange(len(labels)), labels)
    axes[0].invert_yaxis()
    fig.suptitle('Observed profile offsets, independent of the biochemical model', fontsize=12)
    fig.text(.5, .015, 'Positive: intronic crossing occurs earlier. Bars: nominal 95% paired-bootstrap intervals, 5,000 draws.\n'
             'Intervals condition on assigned stages and inferred specimen pairing; no multiplicity correction.',
             ha='center', va='bottom', fontsize=8.5, color='.3')
    fig.tight_layout(rect=(0, .13, 1, .93))
    for ext in ['png', 'svg']:
        fig.savefig(ROOT / f'empirical_profile_offsets.{ext}', dpi=190)
    plt.close(fig)

if __name__ == '__main__':
    main()
