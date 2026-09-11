"""Constrained fits and their residuals; no parameters are estimated here."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pulse_fit import PulseFit, load, GENES, X
from pulse_forward import pulse_values

ROOT = Path(__file__).resolve().parent
COLORS = ['#256a9a', '#ad4937']

def main():
    pairs = load()
    plt.rcParams.update({'font.size': 9, 'svg.fonttype': 'none', 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.linewidth': .7})
    fig, axes = plt.subplots(3, 3, figsize=(13, 9.0), sharex='col',
                              gridspec_kw={'height_ratios': [2.1, 1.2, 1.]})
    rows = []
    curve_rows = []
    dense = np.linspace(1, 24, 1151)
    for col, gene in enumerate(GENES):
        mean = np.array([p.mean(0) for p in pairs[gene]])
        sem = np.array([p.std(0, ddof=1)/np.sqrt(len(p)) for p in pairs[gene]])
        scale = mean.max(0)
        y, se = mean/scale, sem/scale
        fits = {shape: json.loads((ROOT/'fitresults'/f'fit_{gene}_{shape}.json').read_text())
                for shape in ['cosine', 'gaussian']}
        for shape, style in [('cosine', '-'), ('gaussian', '--')]:
            for mode, fit in fits[shape].items():
                pred = np.asarray(fit['pred'])
                for idx in range(len(X)):
                    curve_rows.append(dict(gene=gene, shape=shape, mode=mode,
                        stage_index=int(X[idx]), intron=y[idx,0], exon=y[idx,1],
                        intron_fit=pred[idx,0], exon_fit=pred[idx,1],
                        intron_residual=pred[idx,0]-y[idx,0],
                        exon_residual=pred[idx,1]-y[idx,1]))
                rms = np.sqrt(np.mean((pred-y)**2, axis=0))/np.ptp(y, axis=0)*100
                rows.append(dict(gene=gene, shape=shape, mode=mode, SSE=fit['cost'],
                     intron_RMS_percent_range=rms[0], exon_RMS_percent_range=rms[1],
                     max_parameters=(10 if gene=='Kr' else 14)-(mode=='collapsed'),
                     T=fit['parameters']['T'], processing=fit['parameters']['splice'],
                     decay=fit['parameters']['decay'], gain=fit['parameters']['gain'],
                     bounds='; '.join(fit.get('bound_hits', []))))
            fit = fits[shape]['mapped']
            model = PulseFit(gene, shape=shape, dt=fit['dt'])
            prediction = model.predict(fit, dense)
            for c, color in enumerate(COLORS):
                axes[0,col].plot(dense, prediction[:,c], style, color=color, lw=1.7, alpha=.95)
            initiation = sum(a*pulse_values(dense, q['center'], q['rise'], q['fall'], shape)
                             for a,q in zip(fit['amplitudes'],fit['parameters']['pulses']))
            axes[2,col].plot(dense, initiation/max(initiation.max(),1e-12), style, color='.25', lw=1.5)
        primary = fits['cosine']['mapped']
        for c, color in enumerate(COLORS):
            axes[0,col].errorbar(X, y[:,c], yerr=se[:,c], fmt='o', color=color,
                                  ms=3, lw=.7, capsize=1.5, zorder=5)
            residual = y[:,c]-np.asarray(primary['pred'])[:,c]
            axes[1,col].errorbar(X, residual, yerr=se[:,c], fmt='o-', color=color,
                                  ms=2.8, lw=.7, capsize=1.3)
        axes[0,col].set_title(gene, loc='left', fontstyle='italic', fontweight='bold', fontsize=12)
        axes[0,col].set_ylim(-.08, 1.24)
        axes[1,col].axhline(0, color='.5', ls=':', lw=1)
        axes[2,col].set_ylim(-.03, 1.13)
        for row in range(3):
            axes[row,col].set_xlim(.6,24.4)
            axes[row,col].grid(axis='y', alpha=.1)
        axes[2,col].set_xticks(np.arange(1,25,3), [f'{n}.1' for n in range(1,9)], rotation=45)
        axes[2,col].set_xlabel('eve-en stage (equal substage spacing)')
    axes[0,0].set_ylabel('Signal / observed channel maximum')
    axes[1,0].set_ylabel('Observed minus primary fit')
    axes[2,0].set_ylabel('Initiation / own peak')
    handles = [Line2D([],[],color=COLORS[0],marker='o',lw=0,label='Intron observations'),
               Line2D([],[],color=COLORS[1],marker='o',lw=0,label='Exon observations'),
               Line2D([],[],color='.25',lw=1.6,label='Primary: raised-cosine pulse(s)'),
               Line2D([],[],color='.25',lw=1.6,ls='--',label='Sensitivity: Gaussian pulse(s)')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.966),ncol=4,frameon=False,fontsize=9)
    fig.suptitle('Full mapped-target model with a constrained shared initiation input', fontsize=13, y=.995)
    fig.text(.5,.006,'Both channels are fitted jointly. Error bars: ±1 SEM across 4–5 paired specimens per stage.\n'
             'One initiation pulse for Kr; two for mlpt and svb. Fits are descriptive; individual biochemical rates are not identified.',
             ha='center',va='bottom',fontsize=8.5,color='.3')
    fig.tight_layout(rect=(0,.055,1,.92),h_pad=1.6,w_pad=1.4)
    for ext in ['png','svg']:
        fig.savefig(ROOT/f'constrained_fitted_profiles.{ext}',dpi=190)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(ROOT/'constrained_fit_summary.csv',index=False)
    pd.DataFrame(rows).to_csv(ROOT/'fitresults/fit_summary.csv',index=False)
    pd.DataFrame(curve_rows).to_csv(ROOT/'fitresults/normalized_curves.csv',index=False)

if __name__=='__main__':
    main()
