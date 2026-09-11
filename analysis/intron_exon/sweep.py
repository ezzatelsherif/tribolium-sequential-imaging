"""Prespecified hypothetical-pulse sweep of mapped observation kernels.

See sweep_design.md. Run `python sweep.py` for sweep, summaries, numerical
checks, and figures. No fitting or experimental observations enter this sweep.
"""
from __future__ import annotations
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parent
from position_kernels import cumulative, age_kernel

OUT=ROOT/'sweep_results'
GRID=(.01,.03,.1,.3,1.,3.)
SHAPES=(.25,.5,.75)
TOL=.01
DT=.001


def load_geometry():
    p=ROOT/'data'/'geometry.json'
    return json.loads(p.read_text())['genes']


def pulse(t,r):
    t=np.asarray(t,dtype=float)
    out=np.zeros_like(t)
    left=(t>=0)&(t<=r)
    right=(t>r)&(t<=1)
    out[left]=np.sin(np.pi*t[left]/(2*r))**2
    out[right]=np.cos(np.pi*(t[right]-r)/(2*(1-r)))**2
    return out


def simulate(g,T,S,L,r,dt=DT,mode='mapped'):
    # Kernel age bins and initiation samples both have centers (j+1/2)dt.
    # Their convolution output therefore has physical times (n+1)dt.
    n_age=int(np.ceil((T+12*max(S,L))/dt))
    edges=np.arange(n_age+1)*dt
    f0i,_,f0e,_=cumulative(edges,g,T,S,L,mode)
    wi=np.maximum(np.diff(f0i),0)
    we=np.maximum(np.diff(f0e),0)
    u=pulse((np.arange(int(round(1/dt)))+.5)*dt,r)
    i=np.maximum(fftconvolve(u,wi),0)
    e=np.maximum(fftconvolve(u,we),0)
    t=np.arange(1,len(i)+1)*dt
    # Include known zero at exact time zero to document the boundary.
    t=np.r_[0.,t]; i=np.r_[0.,i]; e=np.r_[0.,e]
    full_i,_,full_e,_=cumulative(np.array([T+200*max(S,L)]),g,T,S,L,mode)
    mass_i=abs(np.trapezoid(i,t)-.5*full_i[0])/(.5*full_i[0])
    mass_e=abs(np.trapezoid(e,t)-.5*full_e[0])/(.5*full_e[0])
    return t,i,e,{'mass_relative_error_I':mass_i,'mass_relative_error_E':mass_e,
                  'tail_fraction_I':i[-1]/i.max(),'tail_fraction_E':e[-1]/e.max()}


def features(t,y):
    y=y/y.max()
    j=int(np.argmax(y))
    peak=t[j]
    if 0<j<len(y)-1:
        den=y[j-1]-2*y[j]+y[j+1]
        if den<0:
            shift=.5*(y[j-1]-y[j+1])/den
            if abs(shift)<=1: peak+=shift*(t[1]-t[0])
    out={'peak':peak,'flat_peak_width':np.ptp(t[y>=1-1e-6])}
    for q in (.4,.5,.6):
        lab=str(int(q*100))
        up=np.where((y[:-1]<q)&(y[1:]>=q))[0]
        down=np.where((y[:-1]>=q)&(y[1:]<q))[0]
        def crossing(k):
            return t[k]+(q-y[k])/(y[k+1]-y[k])*(t[k+1]-t[k])
        out['rise_'+lab]=crossing(up[0]) if len(up) else np.nan
        out['fall_'+lab]=crossing(down[-1]) if len(down) else np.nan
        out['upcrossings_'+lab]=len(up)
        out['downcrossings_'+lab]=len(down)
    return out


def classify(v):
    return 'invalid' if not np.isfinite(v) else 'intron_lead' if v>TOL else 'exon_lead' if v < -TOL else 'small'


def scenario(gene,g,T,S,L,r,mode='mapped',dt=DT):
    t,i,e,diag=simulate(g,T,S,L,r,dt,mode)
    fi,fe=features(t,i),features(t,e)
    row={'gene':gene,'mode':mode,'T':T,'processing_mean':S,'decay_mean':L,
         'peak_fraction':r,'dt':dt,**diag}
    for k,v in fi.items():row['I_'+k]=v
    for k,v in fe.items():row['E_'+k]=v
    row['peak_offset']=fe['peak']-fi['peak']
    for ph,q in itertools.product(('rise','fall'),(40,50,60)):
        key=f'{ph}_{q}'
        row[key+'_offset']=fe[key]-fi[key]
        row[key+'_class']=classify(row[key+'_offset'])
    row['valid']=all(np.isfinite(row[f'{ph}_{q}_offset']) for ph,q in itertools.product(('rise','fall'),(40,50,60)))
    return row


def summarize(df,by):
    rows=[]
    groups=df.groupby(by,dropna=False) if by else [((),df)]
    for key,d in groups:
        if not isinstance(key,tuple):key=(key,)
        for ph,q in itertools.product(('rise','fall'),(40,50,60)):
            vals=d[f'{ph}_{q}_offset']
            classes=d[f'{ph}_{q}_class']
            rec={**dict(zip(by,key)),'phase':ph,'threshold':q,'n_scenarios':len(d),
                 'n_valid':int(vals.notna().sum()),'n_invalid':int(vals.isna().sum()),
                 'n_intron_lead':int((classes=='intron_lead').sum()),
                 'n_exon_lead':int((classes=='exon_lead').sum()),
                 'n_small':int((classes=='small').sum()),
                 'n_positive_without_tolerance':int((vals>0).sum()),
                 'n_negative_without_tolerance':int((vals<0).sum()),
                 'n_zero_without_tolerance':int((vals==0).sum()),
                 'minimum_offset':vals.min(),'maximum_offset':vals.max(),'median_offset':vals.median()}
            rec['fraction_intron_lead']=rec['n_intron_lead']/len(d)
            rec['fraction_exon_lead']=rec['n_exon_lead']/len(d)
            rec['fraction_small']=rec['n_small']/len(d)
            rec['fraction_positive_without_tolerance']=rec['n_positive_without_tolerance']/len(d)
            rec['fraction_negative_without_tolerance']=rec['n_negative_without_tolerance']/len(d)
            rows.append(rec)
    return pd.DataFrame(rows)


def numerical_checks(df,geometries):
    # All 72 mapped corners, plus up to 72 nearest classification/sign boundaries.
    corner=df[df['T'].isin((.01,3.)) & df.processing_mean.isin((.01,3.)) & df.decay_mean.isin((.01,3.))]
    ids=set(corner.index)
    for ph,boundary in itertools.product(('rise','fall'),(-TOL,0.,TOL)):
        ids.update((df[f'{ph}_50_offset']-boundary).abs().nsmallest(12).index)
    rows=[]
    for k,idx in enumerate(sorted(ids)):
        row=df.loc[idx]
        fine=[]
        for dt in (DT/2,DT/4):
            rf=scenario(row.gene,geometries[row.gene],row['T'],row.processing_mean,row.decay_mean,row.peak_fraction,dt=dt)
            rr={'scenario_index':idx,'gene':row.gene,'T':row['T'],'processing_mean':row.processing_mean,
                'decay_mean':row.decay_mean,'peak_fraction':row.peak_fraction,'dt':dt}
            for ph in ('rise','fall'):
                rr[ph+'_50_change']=rf[ph+'_50_offset']-row[ph+'_50_offset']
                rr[ph+'_class_changed']=rf[ph+'_50_class']!=row[ph+'_50_class']
                rr[ph+'_sign_changed']=np.sign(rf[ph+'_50_offset'])!=np.sign(row[ph+'_50_offset'])
            rr['peak_change']=rf['peak_offset']-row['peak_offset']
            tc,ic,ec,_=simulate(geometries[row.gene],row['T'],row.processing_mean,row.decay_mean,row.peak_fraction,DT)
            tf,inf,ef,_=simulate(geometries[row.gene],row['T'],row.processing_mean,row.decay_mean,row.peak_fraction,dt)
            rr['max_normalized_signal_change']=max(np.max(abs(ic/ic.max()-np.interp(tc,tf,inf/inf.max()))),np.max(abs(ec/ec.max()-np.interp(tc,tf,ef/ef.max()))))
            rows.append(rr)
        if (k+1)%30==0: print('Convergence',k+1,'/',len(ids),flush=True)
    check=pd.DataFrame(rows)
    check.to_csv(OUT/'numerical_convergence.csv',index=False)
    # Independent quadrature, three distinct mapped scenarios, near channel
    # maxima and half-height crossings. Split at kernel/pulse breakpoints.
    records=[]
    choices=[('Kr',.01,.03,3.,.25),('mlpt',3.,.01,.1,.5),('svb',.3,3.,.01,.75)]
    for gene,T,S,L,r in choices:
        g=geometries[gene]
        t,i,e,_=simulate(g,T,S,L,r)
        fi,fe=features(t,i),features(t,e)
        for channel,y,f in [('I',i,fi),('E',e,fe)]:
            for feature in ('rise_50','peak','fall_50'):
                obs_t=f[feature]
                breaks=[0.,1.,r,obs_t,obs_t-T,obs_t-T*g['splice_site_bp']/g['transcript_end_bp']]
                for a,b in g['intron_targets']+g['exon_targets']:
                    breaks.extend([obs_t-a*T/g['transcript_end_bp'],obs_t-b*T/g['transcript_end_bp']])
                breaks=sorted(set(v for v in breaks if 0<=v<=1))
                def fun(s):
                    ki,ke=age_kernel(obs_t-s,g,T,S,L)
                    return float(pulse(s,r)*(ki if channel=='I' else ke))
                exact=sum(quad(fun,a,b,epsabs=1e-11,epsrel=1e-10)[0] for a,b in zip(breaks[:-1],breaks[1:]))
                approx=np.interp(obs_t,t,y)
                records.append({'gene':gene,'T':T,'processing_mean':S,'decay_mean':L,'peak_fraction':r,
                    'channel':channel,'feature':feature,'time':obs_t,'quadrature_signal':exact,'binned_signal':approx,
                    'absolute_error':abs(exact-approx),'peak_normalized_error':abs(exact-approx)/y.max()})
    pd.DataFrame(records).to_csv(OUT/'independent_quadrature.csv',index=False)
    return check


def make_figures(df):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(10.3,6.7),sharex=True,sharey=True)
    genes=['Kr','mlpt','svb']
    for j,gene in enumerate(genes):
        for k,ph in enumerate(('rise','fall')):
            d=df[df.gene==gene]
            grid=np.array([[np.mean(d[(d['T']==T)&(d.decay_mean==L)][f'{ph}_50_offset']>TOL) for L in GRID] for T in GRID])
            ax=axes[k,j]
            im=ax.imshow(grid,origin='lower',vmin=0,vmax=1,cmap='Blues',aspect='equal')
            for a,b in itertools.product(range(6),repeat=2):
                n=int(round(grid[a,b]*18))
                ax.text(b,a,str(n),ha='center',va='center',fontsize=8,color='white' if grid[a,b]>.57 else 'black')
            ax.set_xticks(range(6),[str(v) for v in GRID]); ax.set_yticks(range(6),[str(v) for v in GRID])
            ax.set_title(f'{gene}: {"rising" if ph=="rise" else "falling"} half-height')
            if j==0:ax.set_ylabel('Traversal duration / pulse duration')
            if k==1:ax.set_xlabel('Mean RNA decay wait / pulse duration')
    fig.subplots_adjust(left=.085,right=.88,bottom=.12,top=.87,wspace=.15,hspace=.29)
    cb=fig.colorbar(im,cax=fig.add_axes([.9,.2,.017,.54]));cb.set_label('Fraction of the 18 specified scenarios')
    fig.suptitle('When does the intronic profile cross first?',fontsize=15,y=.98)
    fig.text(.5,.925,'Each cell: 6 processing waits x 3 pulse shapes. Numbers show counts out of 18.',ha='center',fontsize=10)
    fig.text(.5,.025,'Exon minus intron crossing > 0.01 pulse duration. Scenario fractions are not biological probabilities.',ha='center',fontsize=9)
    for suffix in ('png','svg'):fig.savefig(ROOT/f'sweep_phase_map.{suffix}',dpi=190)
    plt.close(fig)
    # Signed offsets reveal reversals and near-zero outcomes instead of hiding
    # their magnitudes behind a single percentage.
    fig,axes=plt.subplots(1,3,figsize=(10.4,3.6),sharex=True,sharey=True)
    for ax,gene in zip(axes,genes):
        d=df[df.gene==gene]
        for r,c in zip(SHAPES,('#0072B2','#D55E00','#009E73')):
            s=d[d.peak_fraction==r]
            ax.scatter(s.rise_50_offset,s.fall_50_offset,s=13,c=c,alpha=.48,label=f'Peak at {r:g} D',linewidths=0)
        ax.axvspan(-TOL,TOL,color='grey',alpha=.15);ax.axhspan(-TOL,TOL,color='grey',alpha=.15)
        ax.axvline(0,color='.4',lw=.6);ax.axhline(0,color='.4',lw=.6)
        ax.set_title(gene);ax.set_xlabel('Rising half-height offset / D')
    axes[0].set_ylabel('Falling half-height offset / D')
    axes[-1].legend(frameon=False,fontsize=8,loc='upper left')
    fig.suptitle('All 1,944 mapped-model scenarios',y=.98)
    fig.tight_layout(rect=(0,.08,1,.93))
    fig.text(.5,.015,'D is hypothetical initiation-pulse duration; positive offsets mean the intronic profile crosses first.',ha='center',fontsize=9)
    for suffix in ('png','svg'):fig.savefig(ROOT/f'sweep_signed_offsets.{suffix}',dpi=190)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    geometries=load_geometry()
    rows=[]
    for gene,g in geometries.items():
        for T,S,L,r in itertools.product(GRID,GRID,GRID,SHAPES):
            rows.append(scenario(gene,g,T,S,L,r))
        print('Mapped scenarios complete:',gene,len(rows),flush=True)
    df=pd.DataFrame(rows)
    df.to_csv(OUT/'sweep_results.csv',index=False)
    collapsed=[]
    for S,L,r in itertools.product(GRID,GRID,SHAPES):
        collapsed.append(scenario('geometry_independent',geometries['Kr'],0.,S,L,r,mode='collapsed'))
    dc=pd.DataFrame(collapsed);dc.to_csv(OUT/'collapsed_results.csv',index=False)
    for name,by in [('overall',[]),('by_gene',['gene']),('by_gene_shape',['gene','peak_fraction']),('by_gene_kinetics',['gene','T','processing_mean','decay_mean']),('by_gene_traversal_decay',['gene','T','decay_mean'])]:
        summarize(df,by).to_csv(OUT/f'summary_{name}.csv',index=False)
    summarize(dc,[]).to_csv(OUT/'summary_collapsed.csv',index=False)
    check=numerical_checks(df,geometries)
    make_figures(df)
    print(summarize(df,['gene']).query('threshold==50').to_string(index=False),flush=True)
    print('Max absolute 50% offset changes:',check[['rise_50_change','fall_50_change']].abs().max().to_dict(),flush=True)
    print('Classification changes:',check[['rise_class_changed','fall_class_changed']].sum().to_dict(),flush=True)


if __name__=='__main__':main()
