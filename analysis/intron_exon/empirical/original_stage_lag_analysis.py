"""Exploratory intron/exon analysis. Does not change source workbooks.

Run in this directory with Python (numpy, pandas, scipy, openpyxl, matplotlib).
Pairing within stage follows numbered workbook columns and needs provenance
confirmation. Replicate numbers are NOT treated as longitudinal embryos.
Intervals condition on assigned stages and do not include staging uncertainty.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
SEED = 20260909
N_BOOT = 5000
STAGES = [f'{i}.{j}' for i in range(1, 9) for j in (1, 2, 3)]
X = np.arange(1, 25, dtype=float)
WINDOWS = {'KR': [(1, 12)], 'MLPT': [(1, 13), (13, 24)],
           'SVB': [(1, 12), (12, 24)]}

def load_gene(gene):
    df = pd.read_excel(ROOT / 'temporal' / f'{gene}.xlsx', sheet_name='Sheet1')
    assert [f'{v:.1f}' for v in df.Stage] == STAGES
    # stage, specimen, channel. NaN propagates when either A or B is missing.
    data = np.array([[[r[f'{c}{j}A'] - r[f'{c}{j}B']
                       for c in ('mrna', 'intr')]
                      for j in range(1, 6)] for _, r in df.iterrows()])
    rows = []
    pairs = []
    for s, a in enumerate(data):
        valid = np.isfinite(a).all(axis=1)
        pairs.append(a[valid])
        for j, v in enumerate(a):
            rows.append(dict(gene=gene, stage=STAGES[s], stage_index=s+1,
                             specimen_label=j+1, exon=v[0], intron=v[1],
                             complete_pair=bool(valid[j])))
    return pairs, rows

def crossings(means, window, level=0.5, method='pchip', local_peak=True):
    """Crossings relative to each channel's whole-series minimum and pulse peak.

    This gives both pulses their own half-height; global-peak sensitivity is
    also evaluated because the original display normalized the whole series.
    """
    lo, hi = window
    times = []
    for c in range(2):
        y = means[:, c]
        pk = lo - 1 + int(np.argmax(y[lo-1:hi]))
        top = y[pk] if local_peak else np.max(y)
        threshold = np.min(y) + level * (top - np.min(y))
        z = y - threshold
        rise = [j for j in range(lo-1, pk) if z[j] <= 0 < z[j+1]]
        fall = [j for j in range(pk, hi-1) if z[j] >= 0 > z[j+1]]
        if not rise or not fall:
            times.append([np.nan, np.nan])
            continue
        interp = PchipInterpolator(X, z) if method == 'pchip' else None
        edges = []
        for j in (rise[-1], fall[0]):
            t = (brentq(interp, X[j], X[j+1]) if interp is not None else
                 X[j] - z[j] / (z[j+1] - z[j]))
            edges.append(t)
        times.append(edges)
    onset, decline = np.array(times)[0] - np.array(times)[1]
    return np.array([onset, decline, decline-onset])

def main(output_dir=None, draws=N_BOOT):
    global N_BOOT
    N_BOOT = draws
    OUT = ROOT / 'derived' if output_dir is None else Path(output_dir).resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows, summaries, sensitivity = [], [], []
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8), sharex=True)
    for gi, gene in enumerate(WINDOWS):
        pairs, rows = load_gene(gene)
        all_rows.extend(rows)
        mean = np.array([a.mean(axis=0) for a in pairs])
        rng = np.random.default_rng(SEED + gi)
        # Resample paired specimens independently WITHIN each stage.
        bm = np.stack([a[rng.integers(0, len(a), (N_BOOT, len(a)))].mean(axis=1)
                       for a in pairs], axis=1)
        estimates = np.stack([crossings(mean, w) for w in WINDOWS[gene]])
        bs = np.array([[crossings(m, w) for w in WINDOWS[gene]] for m in bm])
        for wi, w in enumerate(WINDOWS[gene]):
            for k, metric in enumerate(['rising', 'falling', 'falling_minus_rising']):
                b = bs[:, wi, k]
                good = b[np.isfinite(b)]
                q = np.quantile(good, [.025, .975])
                summaries.append(dict(gene=gene, pulse=wi+1, metric=metric,
                                      estimate=estimates[wi,k], lo95=q[0], hi95=q[1],
                                      valid_bootstraps=len(good), total_bootstraps=N_BOOT))
            for method in ['pchip', 'linear']:
                for level in [.4, .5, .6]:
                    for lp in [True, False]:
                        v = crossings(mean, w, level, method, lp)
                        sensitivity.append(dict(gene=gene,pulse=wi+1,method=method,
                            level=level,peak='pulse' if lp else 'whole_series',
                            rising=v[0],falling=v[1],falling_minus_rising=v[2]))
        # Whole-series min-max display, matching notebook convention.
        norm = (mean-mean.min(axis=0))/(mean.max(axis=0)-mean.min(axis=0))
        bn = (bm-bm.min(axis=1,keepdims=True))/(bm.max(axis=1,keepdims=True)-bm.min(axis=1,keepdims=True))
        for c, (color, label) in enumerate([('#b34335','Exon'),('#4275aa','Intron')]):
            dense = np.linspace(1,24,1200)
            axes[gi].plot(dense,PchipInterpolator(X,norm[:,c])(dense),color=color,label=label)
            lower,upper=np.quantile(bn[:,:,c],[.025,.975],axis=0)
            axes[gi].errorbar(X,norm[:,c],yerr=[norm[:,c]-np.minimum(lower,norm[:,c]),
                            np.maximum(upper,norm[:,c])-norm[:,c]],fmt='o',ms=3,
                            color=color,lw=.8,capsize=2)
        axes[gi].set_ylabel(f'{gene}\nRelative intensity')
        axes[gi].set_ylim(-.08,1.15)
        axes[gi].grid(alpha=.17)
        if len(WINDOWS[gene]) > 1:
            axes[gi].axvline(WINDOWS[gene][1][0],color='.6',ls=':',lw=1)
        print(gene,'complete pairs:',sum(map(len,pairs)),flush=True)
    pd.DataFrame(all_rows).to_csv(OUT/'background_corrected_pairs.csv',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'stage_lag_bootstrap.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(OUT/'stage_lag_sensitivity.csv',index=False)
    axes[0].legend(loc='upper right',frameon=False)
    axes[-1].set_xticks(X,STAGES,rotation=60)
    axes[-1].set_xlabel('eve-en stage (one unit = one consecutive substage)')
    fig.suptitle('Exploratory reconstruction from raw A-B measurements\nPoints: stage means; error bars: pointwise 95% paired-bootstrap intervals',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95))
    fig.savefig(OUT/'raw_temporal_profiles.png',dpi=180)
    fig.savefig(OUT/'raw_temporal_profiles.pdf')
    print(pd.DataFrame(summaries).round(3).to_string(index=False))
    print('Sensitivity estimates')
    print(pd.DataFrame(sensitivity).groupby(['gene','pulse'])[['rising','falling','falling_minus_rising']].agg(['min','max']).round(3))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Recompute paired stage-profile offsets from source workbooks.')
    parser.add_argument('--output-dir', type=Path, default=None)
    parser.add_argument('--draws', type=int, default=5000)
    args = parser.parse_args()
    if args.draws < 2:
        parser.error('--draws must be at least 2')
    main(args.output_dir, args.draws)
