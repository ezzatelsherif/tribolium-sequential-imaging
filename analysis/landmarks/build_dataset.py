"""Build signed temporal/spatial distances from the reviewed smoothed measurements.

Replay: python build_dataset.py
First extraction accepts a complete peak-analysis work directory via --source-work.
All regenerated features use the SAME saved paired-channel bootstrap draws. No
embryo bootstrap is rerun. Cached feature arrays include observed draw zero.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, shutil, sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'; DATA.mkdir(exist_ok=True)
STAGES={'Stage3_3':8,'Stage4_1':9,'Stage4_3':11,'Stage5_2':13}
GAP={'kr','gt','mlpt','svb'}
PR={'eve','run','odd'}

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def feature_key(co,event,landmark):return f'{co}__{event}__{landmark}'
def canonical_co(co):return 'pairrule' if co in ['pairrule','shared_eve_run_odd'] else co
def group(a,b):return 'gap-gap' if a in GAP and b in GAP else 'pair-rule' if a in PR and b in PR else 'gap-pair-rule'
def truth(x):return str(x).lower() in ['true','1','1.0']
def summary(v):
    good=v[1:][np.isfinite(v[1:])]
    ci=np.quantile(good,[.025,.975]) if len(good) else [np.nan,np.nan]
    return float(v[0]),float(ci[0]),float(ci[1]),len(good)/len(v[1:])
def label(event):
    return event.replace('_cycle',' cycle ').replace('_first',' first').replace('_peak',' peak ').replace('_episode',' episode')

def save_npz_atomic(path, arrays):
    temporary=path.with_name(path.name+'.tmp')
    with temporary.open('wb') as stream:np.savez_compressed(stream,**arrays)
    with np.load(temporary) as check:
        assert set(check.files)==set(arrays)
        for key in arrays:assert np.array_equal(check[key],arrays[key],equal_nan=True),key
    temporary.replace(path)

def prepare_sources(work):
    needed={
      'peak_inventory.csv':'figure7_ver05/Peak_Order_Comparisons_ver05.csv',
      'posterior_mask.csv':'figure7_ver05_analysis/pairrule/posterior_mask.csv',
      'kr_spatial.csv':'figure7_ver04_analysis/kr_eve/spatial_features.csv',
      'svb_spatial.csv':'figure7_ver04_analysis/svb_eve/spatial_features.csv',
      'gt_mlpt_spatial.csv':'figure7_ver04_analysis/gt_mlpt/spatial_domain_metrics.csv',
      'pairrule_spatial.csv':'figure7_ver04_analysis/pairrule/spatial_features_gene_episode.csv',
      'old_gt_mlpt_temporal_borders.csv':'figure7_ver04_analysis/gt_mlpt/temporal_gap_events.csv',
      'old_kr_eve_temporal.csv':'figure7_ver04_analysis/kr_eve/temporal_features.csv',
      'old_svb_eve_temporal.csv':'figure7_ver04_analysis/svb_eve/temporal_features.csv',
    }
    dest=DATA/'source_tables';dest.mkdir(exist_ok=True)
    for name,rel in needed.items():
        if not (dest/name).exists():shutil.copyfile(work/rel,dest/name)
    prior=dest/'previous_border_inventory.csv'
    prior_source=work.parents[1]/'summary/halfheight_comparison_inventory.csv'
    if not prior.exists() and prior_source.exists():shutil.copyfile(prior_source,prior)
    return {name:pd.read_csv(dest/name) for name in needed}

def features(work,repo):
    cache=DATA/'temporal_features_with_observed.npz'
    if cache.exists():
        z=np.load(cache);return {k:z[k] for k in z.files}
    f={}
    for suite in ['pairrule','gap_peaks']:
        z=np.load(work/f'figure7_ver05_analysis/{suite}/temporal_peak_bootstraps.npz')
        for k in z.files:
            co,ev=('pairrule',k) if suite=='pairrule' else k.split('__',1)
            f[feature_key(co,ev,'peak')]=z[k]
    old=work/'figure7_ver04_analysis'
    for gene in ['kr','svb']:
        z=np.load(old/f'{gene}_eve/temporal_feature_bootstrap.npz')
        for lm in ['rise','fall']:f[feature_key(gene.upper(),gene+'_first',lm)]=z[f'{gene}_{lm}']
    z=np.load(old/'pairrule/temporal_feature_bootstraps.npz')
    for gene in ['eve','run','odd']:
        for cycle in [2,3,4,5]:
            for lm in ['rise','fall']:f[feature_key('pairrule',f'{gene}_cycle{cycle}',lm)]=z[f'{gene}_cycle{cycle}_{lm}']
    engine=repo/'analysis/peak_order/engine/figure7_ver04_analysis'
    sys.path.insert(0,str(engine))
    import smoothing
    pairrule=module(engine/'pairrule/extract_pairrule.py','original_pairrule')
    gap=module(engine/'gt_mlpt/reanalyze.py','original_gap')
    raw=np.load(work/'all_gene_direct_analysis/temporal_bootstrap_curves.npz')
    # Preserve the joint resampling state by reusing paired-channel draws for all genes.
    for co in ['GT','MLPT']:
        arr=raw[f'paired_{co}_{co}']
        vals={lm:np.empty(len(arr)) for lm in ['rise','fall']}
        for i in range(0,len(arr),500):
            x,y=smoothing.temporal(arr[i:i+500]); out=gap.measurements(x,y,co)
            for lm in vals:vals[lm][i:i+len(y)]=out[lm]
        event='gt_episode' if co=='GT' else 'mlpt_first'
        for lm,v in vals.items():f[feature_key(co,event,lm)]=v
        print('Extracted coherent paired borders',co,flush=True)
    # Every target keeps its own matched eve temporal cohort. Apply the same
    # pair-rule cycle/trough definition across these matched cohorts.
    for co in ['KR','SVB','MLPT','GT']:
        arr=raw[f'paired_{co}_eve']
        vals={(cy,lm):np.empty(len(arr)) for cy in [1,2,3,4,5] for lm in ['rise','fall']}
        for i in range(0,len(arr),500):
            x,y=smoothing.temporal(arr[i:i+500])
            for cy in [1,2,3,4,5]:
                ev,_=pairrule.temporal_events(x,y,cy)
                for lm in ['rise','fall']:vals[(cy,lm)][i:i+len(y)]=ev[lm]
        for (cy,lm),v in vals.items():f[feature_key(co,f'eve_cycle{cy}',lm)]=v
        print('Extracted matched eve borders',co,flush=True)
    for k,v in f.items():assert v.shape==(20001,), (k,v.shape)
    save_npz_atomic(cache,f)
    # Endpoint differences relative to historical unpaired GT/MLPT border cache.
    oldz=np.load(old/'gt_mlpt/temporal_bootstrap_features.npz');audit=[]
    for co in ['GT','MLPT']:
        ev='gt_episode' if co=='GT' else 'mlpt_first'
        for lm in ['rise','fall']:
            a=summary(oldz[f'{co}_{lm}']);b=summary(f[feature_key(co,ev,lm)])
            audit.append(dict(cohort=co,landmark=lm,old_estimate=a[0],new_estimate=b[0],old_low=a[1],new_low=b[1],old_high=a[2],new_high=b[2],max_ci_endpoint_change=max(abs(a[1]-b[1]),abs(a[2]-b[2]))))
    pd.DataFrame(audit).to_csv(DATA/'paired_bootstrap_border_audit.csv',index=False)
    return f

def outcome(reasons,lo,hi,ds):
    if reasons:return 'excluded'
    sign=1 if lo>0 else -1 if hi<0 else 0
    if sign==0:return 'unresolved_temporal'
    if any(abs(v)<1e-12 for v in ds):return 'unresolved_spatial'
    if all(np.sign(v)==sign for v in ds):return 'preserved_all3'
    if all(np.sign(v)==-sign for v in ds):return 'reversed_all3'
    return 'mixed'

def audits(d,f,tables):
    audit=[]
    for row in tables['old_gt_mlpt_temporal_borders.csv'].to_dict('records'):
        co,lm=row['gene'],row['event'];ev='gt_episode' if co=='GT' else 'mlpt_first'
        new=summary(f[feature_key(co,ev,lm)])
        audit.append(dict(cohort=co,landmark=lm,old_estimate=row['estimate'],new_estimate=new[0],
            old_low=row['ci95_low'],new_low=new[1],old_high=row['ci95_high'],new_high=new[2],
            max_ci_endpoint_change=max(abs(row['ci95_low']-new[1]),abs(row['ci95_high']-new[2]))))
    pd.DataFrame(audit).to_csv(DATA/'paired_bootstrap_border_audit.csv',index=False)
    audit=[]
    for co in ['KR','SVB']:
        prior=tables[f'old_{co.lower()}_eve_temporal.csv']
        for lm in ['rise','fall']:
            row=prior[(prior.gene=='eve')&(prior.event==lm)].iloc[0];new=summary(f[feature_key(co,'eve_cycle3',lm)])
            audit.append(dict(cohort=co,landmark=lm,old_estimate=row['estimate'],new_estimate=new[0],
                old_recovery=row['valid_bootstrap_fraction'],new_recovery=new[3],estimate_change=new[0]-row['estimate']))
    pd.DataFrame(audit).to_csv(DATA/'matched_eve_definition_audit.csv',index=False)
    prior=DATA/'source_tables/previous_border_inventory.csv'
    if prior.exists():
        audit=[]
        for x in pd.read_csv(prior).to_dict('records'):
            q=d[(d.group==x['group'])&(d.stage==x['stage'])&(d.landmark==x['landmark']) & (((d.gene_a==x['gene_a'])&(d.gene_b==x['gene_b']))|((d.gene_a==x['gene_b'])&(d.gene_b==x['gene_a'])))]
            if x['group']=='pair-rule':q=q[q.event_a.str.endswith('cycle'+str(x['cycle']))]
            assert len(q)==1,(x,len(q))
            y=q.iloc[0];sgn=1 if y.gene_a==x['gene_a'] else -1;ci=sorted([sgn*y.ci_low,sgn*y.ci_high])
            audit.append(dict(comparison_id=y.comparison_id,old_outcome=x['outcome'],new_outcome=y.outcome,
                changed=x['outcome']!=y.outcome,old_dt=x['temporal_b_minus_a'],new_dt=sgn*y['dt'],
                old_ci_low=x['ci95_low'],new_ci_low=ci[0],old_ci_high=x['ci95_high'],new_ci_high=ci[1]))
        pd.DataFrame(audit).to_csv(DATA/'previous_border_comparison_audit.csv',index=False)

def main(work,repo):
    tables=prepare_sources(work);f=features(work,repo)
    inv=tables['peak_inventory.csv'];masks=tables['posterior_mask.csv'].groupby(['stage','embryo'])[['verified_ap_lo','verified_ap_hi']].first()
    space={}
    for gene in ['kr','svb']:
        for r in tables[f'{gene}_spatial.csv'].query('gene==@gene').to_dict('records'):
            for lm in ['rise','fall']:
                space[(r['stage'],int(r['embryo']),gene+'_first',lm)]=dict(position=r[lm],peak=r['peak'],invalid=truth(r['peak_at_window_edge']),flags='')
    for r in tables['gt_mlpt_spatial.csv'].to_dict('records'):
        ev='gt_episode' if r['gene']=='gt' else 'mlpt_first'
        for lm in ['rise','fall']:space[(r['stage'],int(r['embryo']),ev,lm)]=dict(position=r[lm],peak=r['peak'],invalid=False,flags='')
    pr=tables['pairrule_spatial.csv'];pr=pr[pr.boundary_definition=='outer_episode_envelope']
    for r in pr.to_dict('records'):
        if r['event'] not in ['rise','fall']:continue
        invalid=truth(r['common_peak_clipped']) or truth(r['posterior_developing_cycle']) or not truth(r['peak_inside_gene_window'])
        space[(r['stage'],int(r['embryo']),f"{r['gene']}_cycle{int(r['cycle'])}",r['event'])]=dict(position=r['estimate'],peak=r['common_peak_estimate'],invalid=invalid,flags=str(r['flags']))
    records=[];contrasts={};seen_borders=set()
    for keys,z in inv.groupby(['stage','gene_a','gene_b','event_a','event_b','cohort_a','cohort_b'],sort=False):
        stage,a,b,ea,eb,ca,cb=keys;ca=canonical_co(ca);cb=canonical_co(cb);cls=group(a,b)
        if cls=='pair-rule' and ea.split('cycle')[-1]!=eb.split('cycle')[-1]:continue
        z=z.sort_values('embryo')
        for lm in ['peak','rise','fall']:
            fa=ea if lm=='peak' else ('gt_episode' if a=='gt' else ea)
            fb=eb if lm=='peak' else ('gt_episode' if b=='gt' else eb)
            unique=(stage,ca,fa,cb,fb,lm)
            if unique in seen_borders:continue
            seen_borders.add(unique)
            ka,kb=feature_key(ca,fa,lm),feature_key(cb,fb,lm)
            av=f.get(ka,np.full(20001,np.nan));bv=f.get(kb,np.full(20001,np.nan))
            delta=bv-av;dt,lo,hi,recovery=summary(delta)
            ah=summary(av)[2];bh=summary(bv)[2]
            cid=f'{ca}_{fa}__vs__{cb}_{fb}__{lm}'
            contrasts[cid]=delta[1:]
            reasons=[];flags=[];ds=[];positions={};spans=[]
            if not np.isfinite([av[0],bv[0]]).all():reasons.append('missing_temporal_landmark')
            elif not np.isfinite([ah,bh]).all() or max(ah,bh)>=STAGES[stage]-.5:reasons.append('temporal_landmark_not_secured_before_current_substage')
            if lm=='peak':
                if len(z)!=3 or z.embryo.nunique()!=3:reasons.append('not_three_spatial_embryos')
                if not z.eligible.map(truth).all():
                    reasons.extend(str(s) for v in z.exclusion_reasons.fillna('') for s in v.split(';') if s)
                    if not reasons:reasons.append('source_spatial_exclusion')
                flags.extend(str(s) for v in z['flags'].fillna('') for s in v.split(';') if s)
                for emb in [1,2,3]:
                    zz=z[z.embryo==emb]
                    if len(zz):
                        r=zz.iloc[0];x,y=r.position_a,r.position_b
                        span=1 if r.span_a_hi<r.span_b_lo else -1 if r.span_b_hi<r.span_a_lo else 0
                    else:x=y=np.nan;span=0
                    ds.append(y-x);positions[f'position_a_emb{emb}']=x;positions[f'position_b_emb{emb}']=y;spans.append(span)
            else:
                if stage=='Stage5_2' and ({a,b}&{'kr','svb'}):reasons.append('excluded_late_Kr_or_svb_domain')
                for emb in [1,2,3]:
                    masklo,maskhi=masks.loc[(stage,emb)];xy=[]
                    for side,event in [('a',fa),('b',fb)]:
                        sr=space.get((stage,emb,event,lm),dict(position=np.nan,peak=np.nan,invalid=True,flags='missing_source_event'))
                        val=sr['position'];xy.append(val);positions[f'position_{side}_emb{emb}']=val
                        if not np.isfinite(val):reasons.append('missing_spatial_crossing')
                        elif not masklo<=val<=maskhi:reasons.append('spatial_crossing_outside_posterior_region')
                        if sr['invalid'] or not np.isfinite(sr['peak']):reasons.append('spatial_domain_identity_unresolved')
                        elif not masklo<=sr['peak']<=maskhi:reasons.append('spatial_domain_peak_outside_posterior_region')
                        flags.extend(s for s in sr['flags'].split(';') if s and s!='nan')
                    ds.append(xy[1]-xy[0]);spans.append(int(np.sign(xy[1]-xy[0])) if np.isfinite(xy).all() else 0)
            if not np.isfinite(ds).all():reasons.append('missing_spatial_landmark')
            reasons=sorted(set(reasons));flags=sorted(set(flags))
            o=outcome(reasons,lo,hi,ds)
            sign=1 if lo>0 else -1 if hi<0 else 0
            robust=('unresolved_temporal' if sign==0 else 'preserved_all3' if all(s==sign for s in spans) else 'reversed_all3' if all(s==-sign for s in spans) else 'split_ambiguous_or_mixed') if not reasons else 'excluded'
            records.append(dict(comparison_id=stage+'__'+cid,contrast_id=cid,group=cls,landmark=lm,stage=stage,
                gene_a=a,gene_b=b,event_a=fa,event_b=fb,cohort_a=ca,cohort_b=cb,
                pair_label=label(fa)+' / '+label(fb),time_a=av[0],time_b=bv[0],dt=dt,ci_low=lo,ci_high=hi,
                time_a_ci_high=ah,time_b_ci_high=bh,matched_temporal_channels=ca==cb,
                ds_emb1=ds[0],ds_emb2=ds[1],ds_emb3=ds[2],**positions,
                eligible=not reasons,exclusion_reasons=';'.join(reasons),flags=';'.join(flags),recovery=recovery,
                outcome=o,span_outcome=robust,
                orientation=1 if dt>=0 else -1,
                peak_definition='split_eve_midpoint' if lm=='peak' and any('midpoint' in x or 'split' in x for x in flags) else 'single_peak' if lm=='peak' else 'outer_episode_halfheight',
                temporal_source_a=ka,temporal_source_b=kb))
    d=pd.DataFrame(records);assert not d.comparison_id.duplicated().any()
    d.to_csv(DATA/'landmark_distances.csv',index=False)
    save_npz_atomic(DATA/'temporal_contrast_draws.npz',contrasts)
    summary_rows=[]
    for (g,l),zz in d.groupby(['group','landmark']):
        s=dict(group=g,landmark=l,candidates=len(zz),eligible=int(zz.eligible.sum()),excluded=int((~zz.eligible).sum()))
        for o in ['preserved_all3','reversed_all3','mixed','unresolved_temporal','unresolved_spatial']:s[o]=int((zz.outcome==o).sum())
        summary_rows.append(s)
    pd.DataFrame(summary_rows).to_csv(DATA/'coverage_summary.csv',index=False)
    feature_rows=[]
    for k,v in f.items():
        obs,l,h,rec=summary(v);feature_rows.append(dict(source_id=k,estimate=obs,ci_low=l,ci_high=h,recovery=rec))
    pd.DataFrame(feature_rows).to_csv(DATA/'temporal_feature_summary.csv',index=False)
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print('All-three reversals:')
    print(d[d.outcome=='reversed_all3'][['group','landmark','stage','pair_label','dt','ci_low','ci_high','ds_emb1','ds_emb2','ds_emb3']].to_string(index=False))
    manifest={'source_commit':'a09fafea149c097516646d0ad8bec3d5b0ad8169','draws':20000,
      'canonical_sign':'b minus a; AP increases anterior to posterior','orientation':'fixed observed temporal sign; never reoriented per draw',
      'strict_completion':'both marginal 95% CI upper bounds < stage index - 0.5, all landmarks/classes',
      'source_table_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (DATA/'source_tables').glob('*.csv')},
      'comparison_rows':len(d),'unique_contrasts':len(contrasts),'gene_pairs':len(d[['gene_a','gene_b']].apply(lambda r:tuple(sorted(r)),axis=1).unique())}
    (DATA/'dataset_manifest.json').write_text(json.dumps(manifest,indent=2))
    audits(d,f,tables)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-work',type=Path,default=ROOT.parents[1]/'outputs/peak_order/work');p.add_argument('--repository',type=Path,default=ROOT.parents[1]);a=p.parse_args();main(a.source_work,a.repository)
