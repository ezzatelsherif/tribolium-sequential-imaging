"""Class-neutral, peak-only time-to-space screen: target/EVE and cross-cohort pairs.

Keeps original smoothing and original target peak identities. New EVE temporal
references use identical cycle windows across cohorts. No source file is edited.
"""
from pathlib import Path
import sys, json, itertools, hashlib
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / 'figure7_ver04_analysis'
DATA = ROOT / 'all_gene_direct_analysis'
PR = HERE.parent / 'pairrule'
sys.path.insert(0, str(OLD))
import smoothing

LABELS = [f'{c}.{s}' for c in range(1,9) for s in range(1,4)]
TARGETS = {
    'KR': [('kr_first', 0, 11)],
    'SVB': [('svb_first', 0, 11)],
    'MLPT': [('mlpt_first', 2, 12)],
    'GT': [('gt_peak1', 6, 9), ('gt_peak2', 9, 11)],
}
STAGES = {'Stage3_3': 8, 'Stage4_1': 9, 'Stage4_3': 11, 'Stage5_2': 13}

def summary(a):
    a = np.asarray(a); good = np.isfinite(a[1:]); v = a[1:][good]
    ci = np.quantile(v, [.025,.975]) if np.isfinite(a[0]) and len(v) else [np.nan,np.nan]
    return dict(estimate=float(a[0]), ci95_low=float(ci[0]), ci95_high=float(ci[1]),
                valid_bootstrap_n=int(good.sum()), valid_bootstrap_fraction=float(good.mean()))

def peak_in_window(x,y,lo,hi,choose_interior=False):
    ids=np.flatnonzero((x>=lo)&(x<=hi)); r=np.arange(len(y))
    loc=np.zeros_like(y,dtype=bool)
    loc[:,1:-1]=(y[:,1:-1]>=y[:,:-2])&(y[:,1:-1]>=y[:,2:])&((y[:,1:-1]>y[:,:-2])|(y[:,1:-1]>y[:,2:]))
    if choose_interior:
        k=ids[np.argmax(np.where(loc[:,ids],y[:,ids],-np.inf),axis=1)]
        valid=loc[:,ids].any(axis=1)
    else:
        k=ids[np.argmax(y[:,ids],axis=1)];valid=loc[r,k]
    pos=x[k].copy();pos[~valid]=np.nan
    return pos

def temporal():
    original=np.load(DATA/'temporal_bootstrap_curves.npz')
    features={};rows=[];profiles=[];candidates=[];observed=[]
    stages=pd.read_csv(DATA/'temporal_paired_stage_summary.csv')
    for co in TARGETS:
        for gene in [co,'eve']:
            raw=original[f'paired_{co}_{gene}'];n=len(raw)
            specs=TARGETS[co] if gene==co else [(f'eve_cycle{c}',3*(c-1)-.5,3*(c-1)+2.5) for c in range(1,6)]
            arr={event:np.empty(n) for event,_,_ in specs}
            for start in range(0,n,500):
                x,y=smoothing.temporal(raw[start:start+500],normalize=False)
                for event,lo,hi in specs:
                    arr[event][start:start+len(y)]=peak_in_window(x,y,lo,hi,choose_interior=gene=='eve')
                if start==0:
                    yy=y[0]; baseline=yy.min();amplitude=yy.max()-baseline
                    profiles.append(pd.DataFrame(dict(cohort=co,gene=gene.lower(),substage_index=x,smoothed_intensity=yy,
                         normalized_intensity=(yy-baseline)/amplitude,normalization='full_smoothed_series_range')))
                    peaks,prop=find_peaks(yy,prominence=0)
                    for k,prom in zip(peaks,prop['prominences']):
                        assignments=[ev for ev,lo,hi in specs if lo<=x[k]<=hi and abs(arr[ev][0]-x[k])<1e-8]
                        candidates.append(dict(cohort=co,gene=gene.lower(),position=x[k],height=yy[k],
                            prominence_full_range=prom/amplitude,assigned_events=';'.join(assignments),
                            occurs_by_last_spatial_stage=x[k]<=13,
                            flag='assigned' if assignments else 'unassigned_candidate'))
            for event,lo,hi in specs:
                key=f'{co}__{event}';features[key]=arr[event]
                rows.append(dict(cohort=co,gene=gene.lower(),event_id=event,bootstrap_key=key,
                    cycle=int(event.split('cycle')[1]) if 'cycle' in event else np.nan,
                    window_lo=lo,window_hi=hi,**summary(arr[event])))
            ss=stages[(stages.cohort==co)&(stages.gene==gene)].copy()
            ss['gene']=gene.lower();observed.append(ss)
    pd.DataFrame(rows).to_csv(HERE/'temporal_peak_references.csv',index=False)
    pd.concat(profiles).to_csv(HERE/'temporal_profiles.csv',index=False)
    pd.concat(observed).to_csv(HERE/'temporal_observed_stage_means.csv',index=False)
    pd.DataFrame(candidates).to_csv(HERE/'temporal_candidate_inventory.csv',index=False)
    np.savez_compressed(HERE/'temporal_peak_bootstraps.npz',**features)
    # Exact compatibility with established observed v04 target coordinates.
    previous=[]
    for co,folder in [('KR','kr_eve'),('SVB','svb_eve')]:
        d=pd.read_csv(OLD/folder/'temporal_features.csv')
        r=d[(d.gene==co.lower())&(d.event=='peak')].iloc[0]
        new=next(z for z in rows if z['cohort']==co and z['gene']==co.lower())
        previous.append(dict(event_id=new['event_id'],old=float(r.estimate),new=new['estimate']))
    d=pd.read_csv(OLD/'gt_mlpt'/'temporal_gap_peaks.csv')
    for co,old_event,new_event in [('GT','pulse1','gt_peak1'),('GT','pulse2','gt_peak2'),('MLPT','first','mlpt_first')]:
        r=d[(d.gene==co)&(d.event==old_event)].iloc[0]
        new=next(z for z in rows if z['cohort']==co and z['event_id']==new_event)
        previous.append(dict(event_id=new_event,old=float(r.estimate),new=new['estimate']))
    assert all(abs(z['old']-z['new'])<1e-10 for z in previous)
    (HERE/'temporal_compatibility_check.json').write_text(json.dumps(previous,indent=2))

def spatial_targets():
    rows=[];profiles=[];candidates=[]
    for gene,folder in [('kr','kr_eve'),('svb','svb_eve')]:
        d=pd.read_csv(OLD/folder/'spatial_features.csv')
        for r in d[d.gene==gene].itertuples():
            rows.append(dict(stage=r.stage,embryo=int(r.embryo),gene=gene,event_id=f'{gene}_first',
                cycle=np.nan,position=r.peak,lobe_lo=r.peak,lobe_hi=r.peak,principal_peak=r.peak,
                reference_kind='single_peak',spatial_resolved=bool(np.isfinite(r.peak) and not r.peak_at_window_edge),
                excluded_neural=bool(r.excluded_neural),identity_status='original_first_episode',
                window_lo=r.window_lo,window_hi=r.window_hi,flags='neural_expression' if r.excluded_neural else ''))
        pp=pd.read_csv(OLD/folder/'spatial_profiles.csv');pp=pp[pp.gene==gene].copy()
        profiles.append(pp[['stage','embryo','gene','x_ap','raw_intensity','smoothed_intensity']])
    d=pd.read_csv(OLD/'gt_mlpt'/'spatial_episode_metrics.csv')
    for r in d.itertuples():
        event=r.event if r.gene=='gt' else 'mlpt_first'
        rows.append(dict(stage=r.stage,embryo=int(r.embryo),gene=r.gene,event_id=event,
            cycle=np.nan,position=r.peak,lobe_lo=r.peak,lobe_hi=r.peak,principal_peak=r.peak,
            reference_kind='single_peak',spatial_resolved=bool(r.is_local_maximum and not r.peak_at_window_edge),
            excluded_neural=False,identity_status='original_assigned_gt_lobe' if r.gene=='gt' else 'original_first_episode',
            window_lo=r.window_lo,window_hi=r.window_hi,flags=''))
    pp=pd.read_csv(OLD/'gt_mlpt'/'spatial_profiles.csv')
    pp['smoothed_intensity']=pp.smoothed_original_scale
    profiles.append(pp[['stage','embryo','gene','x_ap','raw_intensity','smoothed_intensity']])
    allp=pd.concat(profiles,ignore_index=True)
    norm=[]
    for (s,e,g),p in allp.groupby(['stage','embryo','gene']):
        p=p.sort_values('x_ap').copy();x=p.x_ap.to_numpy();y=p.smoothed_intensity.to_numpy();span=np.ptp(y)
        p['normalized_intensity']=(y-y.min())/span;norm.append(p)
        peaks,props=find_peaks(y,prominence=0)
        rr=[r for r in rows if r['stage']==s and r['embryo']==e and r['gene']==g]
        for k,prom in zip(peaks,props['prominences']):
            assigned=[r['event_id'] for r in rr if abs(r['position']-x[k])<1e-8]
            candidates.append(dict(stage=s,embryo=int(e),gene=g,position=x[k],height=y[k],
                prominence_full_range=prom/span,assigned_events=';'.join(assigned),flag='assigned' if assigned else 'unassigned_candidate'))
    pd.DataFrame(rows).to_csv(HERE/'spatial_target_peaks_unmasked.csv',index=False)
    pd.concat(norm).to_csv(HERE/'spatial_profiles.csv',index=False)
    pd.DataFrame(candidates).to_csv(HERE/'spatial_candidate_inventory.csv',index=False)

if __name__=='__main__':
    temporal();spatial_targets()
    print(pd.read_csv(HERE/'temporal_peak_references.csv')[['cohort','event_id','estimate','ci95_low','ci95_high','valid_bootstrap_fraction']].to_string(index=False))
