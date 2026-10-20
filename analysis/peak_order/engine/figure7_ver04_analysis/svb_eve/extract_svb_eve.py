"""Repeat v03 matched SVB/EVE analysis after original smoothing; also KR/EVE3 audit.
Temporal windows and the half-height definition are unchanged from v03.
The shared original-smoothing helper owns smoothing parameters/provenance.
"""
from pathlib import Path
import sys,json,importlib.util,argparse
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
sys.path.insert(0,str(HERE.parent))
import smoothing
import feature_helpers as v03
DATA=BASE/'all_gene_direct_analysis'
LABELS=v03.LABELS

def extract_temporal(x,y,pl,ph,bl,bh):
    """Find dense-curve maxima, exact-window clipping flags and outer crossings."""
    xx=np.unique(np.r_[x,pl,ph,bl,bh]); yy=PchipInterpolator(x,y,axis=1)(xx)
    r=np.arange(len(y)); baseline=np.min(y,axis=1)
    inds=np.flatnonzero((xx>=pl)&(xx<=ph)); k=inds[np.argmax(yy[:,inds],axis=1)]
    p=yy[r,k];level=baseline+.5*(p-baseline);peak=xx[k]
    edge=(k==inds[0])|(k==inds[-1]);local=(k>0)&(k<len(xx)-1)
    local &= (p>=yy[r,np.maximum(0,k-1)])&(p>=yy[r,np.minimum(len(xx)-1,k+1)])
    above=yy>=level[:,None]; seg=xx[:-1][None,:]
    out={'peak':peak.copy()}
    flags={'peak_at_window_edge':edge,'peak_resolved_local_maximum':local}
    for event in ['rise','fall']:
        if event=='rise':
            can=(~above[:,:-1])&above[:,1:]&(seg>=bl)&(seg<peak[:,None])
            j=can.argmax(axis=1);censor=above[:,np.flatnonzero(xx==bl)[0]]
        else:
            can=above[:,:-1]&(~above[:,1:])&(seg>=peak[:,None])&(seg<bh)
            j=len(xx)-2-can[:,::-1].argmax(axis=1);censor=above[:,np.flatnonzero(xx==bh)[0]]
        with np.errstate(invalid='ignore',divide='ignore'):
            val=xx[j]+(level-yy[r,j])/(yy[r,j+1]-yy[r,j])*(xx[j+1]-xx[j])
        val[(~can.any(axis=1))|censor|(p<=baseline)]=np.nan
        out[event]=val;flags[event+'_censored']=~np.isfinite(val)
    for vals in out.values(): vals[~local]=np.nan
    return out,baseline,p,level,flags

def brackets(x):
    if not np.isfinite(x):return ('','')
    return (LABELS[max(0,min(23,int(np.floor(x))))],LABELS[max(0,min(23,int(np.ceil(x))))])

def temporal(cohort):
    outdir=HERE if cohort=='SVB' else HERE.parent/'kr_eve';outdir.mkdir(exist_ok=True)
    curves=np.load(DATA/'temporal_bootstrap_curves.npz');stages=pd.read_csv(DATA/'temporal_paired_stage_summary.csv')
    records=[];profiles=[];observed=[];boot={};flagrows=[]
    for gene,key,pl,ph,bl,bh in [(cohort.lower(),cohort,0,11,0,11),('eve','eve',6,8,5,8)]:
        raw=curves[f'paired_{cohort}_{key}'];B=len(raw)
        features={k:np.empty(B) for k in ['peak','rise','fall']};flags={};normalized=np.empty((B,800),dtype=np.float32)
        baseline=np.empty(B);pk=np.empty(B);level=np.empty(B)
        for start in range(0,B,500):
            end=min(start+500,B)
            x,y=smoothing.temporal(raw[start:end],normalize=False)
            f,base,p,l,ff=extract_temporal(x,y,pl,ph,bl,bh)
            for k,val in f.items():features[k][start:end]=val
            for k,val in ff.items():
                if k not in flags:flags[k]=np.empty(B,dtype=bool)
                flags[k][start:end]=val
            baseline[start:end]=base;pk[start:end]=p;level[start:end]=l
            normalized[start:end]=(y-base[:,None])/(p-base)[:,None]
            if start==0:observed_y=y[0].copy()
        meta=dict(cohort=cohort,gene=gene,episode='first' if key==cohort else 'cycle3',
            peak_window_start=LABELS[pl],peak_window_end=LABELS[ph],boundary_window_start=LABELS[bl],boundary_window_end=LABELS[bh],
            normalization='full_smoothed_series_min_to_smoothed_episode_max',baseline_raw=baseline[0],
            peak_intensity_raw=pk[0],half_height_raw=level[0],threshold_fraction=.5,
            smoother='PCHIP800_then_SG81_poly3',smoothing_source='original_filter_shared_method_extension_to_SVB_matched_cohort' if cohort=='SVB' else 'original_filter_reapplied_to_KR_matched_cohort',
            observed_peak_stage='between '+' and '.join(brackets(features['peak'][0])),
            peak_at_window_edge=bool(flags['peak_at_window_edge'][0]),peak_resolved_local_maximum=bool(flags['peak_resolved_local_maximum'][0]))
        for event,vals in features.items():
            row={**meta,'event':event,'unit':'ordinal_substage',**v03.summarize(vals),'observed_boundary_censored':bool(not np.isfinite(vals[0]))}
            row['bracket_start'],row['bracket_end']=brackets(vals[0]);records.append(row);boot[f'{gene}_{event}']=vals
        for name,vals in flags.items():
            boot[f'{gene}_{name}']=vals
            flagrows.append(dict(cohort=cohort,gene=gene,flag=name,observed=bool(vals[0]),bootstrap_fraction=float(vals[1:].mean())))
        q=np.quantile(normalized[1:],[.025,.975],axis=0)
        profiles.append(pd.DataFrame(dict(cohort=cohort,gene=gene,substage_index=x,mean=observed_y,normalized_mean=normalized[0],
            normalized_ci95_low=q[0],normalized_ci95_high=q[1],baseline_raw=baseline[0],peak_intensity_raw=pk[0],
            in_peak_window=(x>=pl)&(x<=ph),in_boundary_window=(x>=bl)&(x<=bh))))
        ss=stages[(stages.cohort==cohort)&(stages.gene==key)].sort_values('substage_index').copy()
        assert np.allclose(ss['mean'],raw[0]);ss['gene']=gene;ss['normalized_mean']=(ss['mean']-baseline[0])/(pk[0]-baseline[0]);ss['normalized_sem']=ss['sem']/(pk[0]-baseline[0]);ss['in_peak_window']=ss.substage_index.between(pl,ph);ss['in_boundary_window']=ss.substage_index.between(bl,bh)
        observed.append(ss)
    contrast=[]
    for event in ['rise','peak','fall']:
        vals=boot[f'{cohort.lower()}_{event}']-boot[f'eve_{event}']
        contrast.append(dict(contrast=cohort.lower()+'_minus_eve3',event=event,unit='ordinal_substage',**v03.summarize(vals),probability_less_than_zero=float(np.mean(vals[1:]<0)),probability_greater_than_zero=float(np.mean(vals[1:]>0))))
    pd.DataFrame(records).to_csv(outdir/'temporal_features.csv',index=False);pd.concat(profiles).to_csv(outdir/'temporal_profiles.csv',index=False);pd.concat(observed).to_csv(outdir/'temporal_observed_profiles.csv',index=False)
    pd.DataFrame(contrast).to_csv(outdir/'temporal_paired_contrasts.csv',index=False);pd.DataFrame(flagrows).to_csv(outdir/'temporal_feature_flags.csv',index=False)
    np.savez_compressed(outdir/'temporal_feature_bootstrap.npz',**boot)
    if cohort=='SVB' and (BASE/'figure7_ver03_analysis/svb_eve/temporal_features.csv').exists():
        prior=pd.read_csv(BASE/'figure7_ver03_analysis/svb_eve/temporal_features.csv')
        comparison=pd.DataFrame(records).merge(prior[['gene','event','estimate','ci95_low','ci95_high']],on=['gene','event'],suffixes=('_smoothed','_raw'))
        comparison['change_ordinal_substages']=comparison.estimate_smoothed-comparison.estimate_raw
        comparison.to_csv(outdir/'temporal_raw_vs_smooth.csv',index=False)
    print(cohort,pd.DataFrame(contrast).to_string(index=False))

def spatial(cohort):
    outdir=HERE if cohort=='SVB' else HERE.parent/'kr_eve';outdir.mkdir(exist_ok=True)
    raw=pd.read_csv(DATA/'spatial_raw_profiles.csv');anchors=pd.read_csv(DATA/'spatial_eve_anchors.csv')
    rows=[];profiles=[];peaks=[];crosses=[]
    for (stage,emb),part in raw.groupby(['stage','embryo']):
        if emb not in [1,2,3]:continue
        for gene in [cohort.lower(),'eve']:
            if gene==cohort.lower():lo,hi=.2,1.
            else:
                a=anchors[(anchors.stage==stage)&(anchors.embryo==emb)&(anchors.cycle==3)].iloc[0];lo,hi=float(a.window_lo),float(a.window_hi)
            p=part[part.gene==gene].sort_values('x_ap').copy();x=p.x_ap.to_numpy();y=p.raw_intensity.to_numpy()
            if gene=='svb':
                kw={'frac':.10,'polyorder':3};xx=x;ys=smoothing.smooth_signal(y,**kw)
                provenance='documented_function_default_extension_SVB_setting_not_recovered'
            else:
                kw=smoothing.spatial_parameters(stage,emb,gene)
                xx,ys=smoothing.spatial(x,y,stage,emb,gene,normalize=False)
                provenance='exact_original_WDL_gene_stage_embryo_setting'
            vals,pp,cc=v03.spatial_feature(xx,ys,lo,hi)
            native_kw=dict(smoothing_fraction=kw['frac'],smoothing_polyorder=kw['polyorder'],smoothing_window_n=smoothing.window_length(len(y),**kw),smoothing_native_n=len(y),smoothing_provenance=provenance)
            meta=dict(stage=stage,embryo=int(emb),gene=gene,episode='first' if gene==cohort.lower() else 'cycle3',normalization='full_smoothed_profile_min_to_smoothed_episode_max',excluded_neural=(stage=='Stage5_2' and gene=='svb'),**native_kw)
            rows.append({**meta,**vals})
            for z in pp:peaks.append({**meta,**z})
            for z in cc:crosses.append({**meta,**z})
            p['smoothed_intensity']=ys;p['normalized_episode_intensity']=(ys-vals['baseline_raw'])/(vals['peak_intensity_raw']-vals['baseline_raw']);p['normalized_raw_episode_intensity']=(y-vals['baseline_raw'])/(vals['peak_intensity_raw']-vals['baseline_raw']);p['in_window']=(x>=lo)&(x<=hi);p['excluded_neural']=meta['excluded_neural'];profiles.append(p)
    sr=pd.DataFrame(rows);contrasts=[]
    for (stage,emb),part in sr.groupby(['stage','embryo']):
        a=part.set_index('gene')
        for event in ['rise','peak','fall']:
            g=cohort.lower();contrasts.append(dict(stage=stage,embryo=emb,event=event,**{g:float(a.loc[g,event]),'eve':float(a.loc['eve',event]),g+'_minus_eve':float(a.loc[g,event]-a.loc['eve',event])},excluded_neural=cohort=='SVB' and stage=='Stage5_2'))
    sr.to_csv(outdir/'spatial_features.csv',index=False);pd.concat(profiles).to_csv(outdir/'spatial_profiles.csv',index=False);pd.DataFrame(peaks).to_csv(outdir/'spatial_prominent_peaks.csv',index=False);pd.DataFrame(crosses).to_csv(outdir/'spatial_all_halfheight_crossings.csv',index=False);pd.DataFrame(contrasts).to_csv(outdir/'spatial_paired_offsets.csv',index=False)
    if cohort=='SVB' and (BASE/'figure7_ver03_analysis/svb_eve/spatial_features.csv').exists():
        prior=pd.read_csv(BASE/'figure7_ver03_analysis/svb_eve/spatial_features.csv')
        comp=sr.merge(prior[['stage','embryo','gene','peak','rise','fall']],on=['stage','embryo','gene'],suffixes=('_smoothed','_raw'))
        for event in ['rise','peak','fall']:comp[event+'_change_ap']=comp[event+'_smoothed']-comp[event+'_raw']
        comp.to_csv(outdir/'spatial_raw_vs_smooth.csv',index=False)
    print(cohort,sr[['stage','embryo','gene','rise','peak','fall','peak_at_window_edge']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--part',choices=['temporal','spatial','all'],default='all');p.add_argument('--cohort',choices=['SVB','KR','all'],default='all');args=p.parse_args()
    for cohort in (['SVB','KR'] if args.cohort=='all' else [args.cohort]):
        if args.part in ['temporal','all']:temporal(cohort)
        if args.part in ['spatial','all']:spatial(cohort)
