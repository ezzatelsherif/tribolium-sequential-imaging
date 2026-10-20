"""Join independent peak extraction to the posterior-family mask and screen all pairs.

Retains an exclusion row for every planned stage/embryo comparison. Bootstrap
intervals are pointwise exploratory intervals conditional on detected features.
No global ordering is assembled across independently sampled temporal cohorts.
"""
from pathlib import Path
import itertools,json,sys
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
PR=HERE.parent/'pairrule'
from extract_gap_peaks import STAGES,TARGETS,summary

def truth(v):
    return str(v).lower() in ['true','1','1.0']

def value_or(row,key,fallback):
    value=row.get(key,np.nan)
    return float(value) if np.isfinite(value) else float(fallback)

def mask_targets():
    masks=pd.read_csv(PR/'posterior_mask.csv')
    targets=pd.read_csv(HERE/'spatial_target_peaks_unmasked.csv')
    out=[]
    for r in targets.to_dict('records'):
        m=masks[(masks.stage==r['stage'])&(masks.embryo==r['embryo'])]
        mm=m.iloc[0]
        complete=truth(mm.get('mask_complete',True))
        lo=float(mm.get('mask_ap_lo',np.nan))
        if not np.isfinite(lo):lo=float(mm['verified_ap_lo'])
        hi=float(mm.get('mask_ap_hi',1.))
        inside=np.isfinite(r['position']) and lo<=r['position']<=hi
        rank=np.nan
        for family in m.to_dict('records'):
            if np.isfinite(family.get('window_lo',np.nan)) and family['window_lo']<=r['position']<=family['window_hi']:
                rank=family['posterior_rank'];break
        r.update(in_posterior_mask=inside,posterior_rank=rank,mask_complete=complete,
                 verified_mask_lo=lo,mask_hi=hi,posterior_developing=False,
                 mask_status='inside_verified_posterior_region' if inside else ('outside_requested_posterior_region' if complete else 'outside_verified_coverage_third_family_unresolved'))
        out.append(r)
    df=pd.DataFrame(out);df.to_csv(HERE/'spatial_peak_references.csv',index=False)
    candidates=pd.read_csv(HERE/'spatial_candidate_inventory.csv')
    for k,r in candidates.iterrows():
        m=df[(df.stage==r.stage)&(df.embryo==r.embryo)].iloc[0]
        candidates.loc[k,'in_verified_posterior_region']=bool(m.verified_mask_lo<=r.position<=m.mask_hi)
        candidates.loc[k,'mask_complete']=bool(m.mask_complete)
    candidates.to_csv(HERE/'spatial_candidate_inventory.csv',index=False)
    return df

def temporal_status(a,stage):
    v=float(a['estimate']);end=STAGES[stage]
    if not np.isfinite(v):return 'unresolved_peak'
    if v<end-.5:return 'before_current_substage'
    if v<=end+.5:return 'within_current_substage'
    return 'later_substage'

def main():
    spatial_targets=mask_targets()
    spatial_pr=pd.read_csv(PR/'spatial_peak_references.csv')
    spatial=pd.concat([spatial_targets,spatial_pr],ignore_index=True)
    t=pd.read_csv(HERE/'temporal_peak_references.csv')
    pt=pd.read_csv(PR/'temporal_peak_references.csv')
    boots=np.load(HERE/'temporal_peak_bootstraps.npz');pboots=np.load(PR/'temporal_peak_bootstraps.npz')
    look={(r.cohort,r.event_id):r._asdict() for r in t.itertuples(index=False)}
    for r in pt.itertuples(index=False):
        z=r._asdict();z['cohort']='pairrule';z['bootstrap_key']=r.event_id
        z.setdefault('window_lo',3*(int(r.cycle)-1)-.5)
        z.setdefault('window_hi',3*(int(r.cycle)-1)+2.5)
        look[('pairrule',r.event_id)]=z
    targetevents=[(co,event) for co,specs in TARGETS.items() for event,_,_ in specs]
    pairs=[]
    for a in targetevents:
        for cycle in range(1,6):pairs.append((a,(a[0],f'eve_cycle{cycle}'),'matched_target_eve'))
    for a,b in itertools.combinations(targetevents,2):
        if a[0]!=b[0]:pairs.append((a,b,'independent_target_target'))
    for a in targetevents:
        for gene in ['run','odd']:
            for cycle in sorted(pt.cycle.unique()):
                event=f'{gene}_cycle{int(cycle)}'
                if ('pairrule',event) in look:pairs.append((a,('pairrule',event),'independent_target_run_odd'))
    contrasts=[];results=[]
    for a,b,kind in pairs:
        ta,tb=look[a],look[b]
        va=boots[ta['bootstrap_key']]
        vb=pboots[tb['bootstrap_key']] if b[0]=='pairrule' else boots[tb['bootstrap_key']]
        assert len(va)==len(vb)==20001
        delta=vb-va;st=summary(delta)
        direction=1 if st['ci95_low']>0 else (-1 if st['ci95_high']<0 else 0)
        resolved=direction!=0 and np.isfinite(st['estimate'])
        pair_id=f'{a[0]}__{a[1]}__vs__{b[0]}__{b[1]}'
        overlap=not(ta['window_hi']<tb['window_lo'] or tb['window_hi']<ta['window_lo'])
        meta=dict(pair_id=pair_id,event_a=a[1],event_b=b[1],gene_a=ta['gene'],gene_b=tb['gene'],
            cohort_a=a[0],cohort_b=b[0],comparison_kind=kind,matched_temporal_channels=kind=='matched_target_eve',
            temporal_a=ta['estimate'],temporal_b=tb['estimate'],temporal_difference_b_minus_a=st['estimate'],
            temporal_ci95_low=st['ci95_low'],temporal_ci95_high=st['ci95_high'],
            valid_bootstrap_fraction=st['valid_bootstrap_fraction'],valid_bootstrap_n=st['valid_bootstrap_n'],
            temporal_order_sign=direction,temporal_order='a_before_b' if direction==1 else ('b_before_a' if direction==-1 else 'unresolved'),
            overlapping_temporal_search_windows=overlap,
            high_feature_validity=st['valid_bootstrap_fraction']>=.95,
            intervals='pointwise_95_percent_conditional_on_valid_features_not_multiplicity_adjusted')
        contrasts.append(meta)
        for stage,end in STAGES.items():
            for embryo in [1,2,3]:
                pa=spatial[(spatial.stage==stage)&(spatial.embryo==embryo)&(spatial.event_id==a[1])]
                pb=spatial[(spatial.stage==stage)&(spatial.embryo==embryo)&(spatial.event_id==b[1])]
                flags=[];exclusions=[]
                sa=temporal_status(ta,stage);sb=temporal_status(tb,stage)
                if sa=='later_substage' or sb=='later_substage':exclusions.append('temporal_peak_in_later_substage')
                if sa=='within_current_substage' or sb=='within_current_substage':flags.append('temporal_peak_within_current_substage')
                if sa=='unresolved_peak' or sb=='unresolved_peak':exclusions.append('temporal_peak_unresolved')
                if not overlap:flags.append('disjoint_temporal_search_windows')
                if not meta['high_feature_validity']:flags.append('conditional_interval_with_lower_peak_validity')
                if kind!='matched_target_eve':flags.append('independent_temporal_cohorts_conditional_on_common_staging')
                r={**meta,'stage':stage,'embryo':embryo,'temporal_status_a':sa,'temporal_status_b':sb,
                   'temporal_both_before_current_substage':sa==sb=='before_current_substage',
                   'temporal_both_ci_before_current_substage':ta['ci95_high']<end-.5 and tb['ci95_high']<end-.5}
                if not len(pa) or not len(pb):
                    exclusions.append('spatial_peak_identity_missing')
                    r.update(position_a=np.nan,position_b=np.nan,spatial_difference_b_minus_a=np.nan,
                             spatial_order_sign=0,spatial_span_order_sign=0,spatial_reference_robust=False,
                             centre_based_outcome='excluded',robust_outcome='excluded')
                else:
                    aa=pa.iloc[0].to_dict();bb=pb.iloc[0].to_dict()
                    for side,z in [('a',aa),('b',bb)]:
                        if not truth(z.get('spatial_resolved',z.get('valid',False))):exclusions.append('spatial_peak_unresolved_'+side)
                        if not truth(z.get('in_posterior_mask',True)):exclusions.append('outside_verified_posterior_region_'+side)
                        if 'in_verified_spatial_mask' in z and not pd.isna(z['in_verified_spatial_mask']) and not truth(z['in_verified_spatial_mask']):exclusions.append('outside_verified_spatial_mask_'+side)
                        if truth(z.get('posterior_developing',False)):exclusions.append('posterior_developing_domain_'+side)
                        if truth(z.get('excluded_neural',False)):exclusions.append('neural_expression_'+side)
                        if not truth(z.get('mask_complete',True)):flags.append('third_posterior_family_unresolved_but_verified_subset_used')
                    xa=float(aa['position']);xb=float(bb['position'])
                    alo=float(aa.get('lobe_lo',aa.get('span_lo',xa)));ahi=float(aa.get('lobe_hi',aa.get('span_hi',xa)))
                    blo=float(bb.get('lobe_lo',bb.get('span_lo',xb)));bhi=float(bb.get('lobe_hi',bb.get('span_hi',xb)))
                    dd=xb-xa;sgn=int(np.sign(dd)) if np.isfinite(dd) else 0
                    span_sign=1 if ahi<blo else (-1 if bhi<alo else 0)
                    if not span_sign and sgn:flags.append('order_depends_on_split_domain_centre')
                    outcome='unresolved_temporal_order' if not resolved else ('preserved' if sgn==direction else ('reversed' if sgn else 'unresolved_spatial_order'))
                    robust='unresolved_temporal_order' if not resolved else ('preserved' if span_sign==direction else ('reversed' if span_sign else 'within_split_domain_or_tied'))
                    a15=value_or(aa,'fifteenpct_position',xa);b15=value_or(bb,'fifteenpct_position',xb)
                    alo15=value_or(aa,'fifteenpct_lobe_lo',alo);ahi15=value_or(aa,'fifteenpct_lobe_hi',ahi)
                    blo15=value_or(bb,'fifteenpct_lobe_lo',blo);bhi15=value_or(bb,'fifteenpct_lobe_hi',bhi)
                    sign15=int(np.sign(b15-a15)) if np.isfinite(b15-a15) else 0
                    span15=1 if ahi15<blo15 else (-1 if bhi15<alo15 else 0)
                    outcome15='unresolved_temporal_order' if not resolved else ('preserved' if sign15==direction else ('reversed' if sign15 else 'unresolved_spatial_order'))
                    robust15='unresolved_temporal_order' if not resolved else ('preserved' if span15==direction else ('reversed' if span15 else 'within_split_domain_or_tied'))
                    changed=outcome15!=outcome or robust15!=robust
                    if changed:flags.append('comparison_sensitive_to_eve_split_prominence_criterion')
                    r.update(position_a=xa,position_b=xb,lobe_lo_a=alo,lobe_hi_a=ahi,lobe_lo_b=blo,lobe_hi_b=bhi,
                             posterior_rank_a=aa.get('posterior_rank',np.nan),posterior_rank_b=bb.get('posterior_rank',np.nan),
                             reference_kind_a=aa.get('reference_kind','single_peak'),reference_kind_b=bb.get('reference_kind','single_peak'),
                             spatial_difference_b_minus_a=dd,spatial_order_sign=sgn,spatial_span_order_sign=span_sign,
                             spatial_reference_robust=span_sign!=0,centre_based_outcome=outcome,robust_outcome=robust,
                             fifteenpct_centre_outcome=outcome15,fifteenpct_robust_outcome=robust15,
                             split_criterion_changes_centre_outcome=outcome15!=outcome,
                             split_criterion_changes_robust_outcome=robust15!=robust)
                r['eligible']=not exclusions
                r['exclusion_reasons']=';'.join(sorted(set(exclusions)))
                r['flags']=';'.join(sorted(set(flags)))
                if exclusions:
                    r['centre_based_outcome']='excluded';r['robust_outcome']='excluded'
                    r['fifteenpct_centre_outcome']='excluded';r['fifteenpct_robust_outcome']='excluded'
                r['strong_claim_eligible']=bool(r['eligible'] and r['matched_temporal_channels'] and r['high_feature_validity']
                    and r['temporal_both_ci_before_current_substage'] and r['robust_outcome'] in ['preserved','reversed']
                    and r.get('fifteenpct_robust_outcome')==r['robust_outcome'])
                results.append(r)
    cc=pd.DataFrame(contrasts);rr=pd.DataFrame(results)
    cc.to_csv(HERE/'temporal_pairwise_contrasts.csv',index=False)
    rr.to_csv(HERE/'all_pairwise_orders.csv',index=False)
    rr[rr.eligible].to_csv(HERE/'eligible_pairwise_orders.csv',index=False)
    ss=[]
    for (pair,stage),d in rr.groupby(['pair_id','stage'],sort=False):
        first=d.iloc[0];keep=d[d.eligible]
        z={k:first[k] for k in ['pair_id','event_a','event_b','gene_a','gene_b','cohort_a','cohort_b','comparison_kind','matched_temporal_channels','temporal_a','temporal_b','temporal_ci95_low','temporal_ci95_high','valid_bootstrap_fraction','overlapping_temporal_search_windows']}
        z.update(stage=stage,n_observed_embryos=3,n_eligible_embryos=len(keep))
        for label in ['preserved','reversed','unresolved_temporal_order','unresolved_spatial_order']:
            z['n_centre_'+label]=int((keep.centre_based_outcome==label).sum())
        for label in ['preserved','reversed','unresolved_temporal_order','within_split_domain_or_tied']:
            z['n_robust_'+label]=int((keep.robust_outcome==label).sum())
        z['n_temporal_both_before_current_substage']=int(keep.temporal_both_before_current_substage.sum())
        z['n_split_criterion_changes_centre_outcome']=int((keep.split_criterion_changes_centre_outcome==True).sum())
        z['n_split_criterion_changes_robust_outcome']=int((keep.split_criterion_changes_robust_outcome==True).sum())
        z['n_strong_claim_eligible']=int(keep.strong_claim_eligible.sum())
        z['all_three_robust_outcome']=('preserved' if z['n_robust_preserved']==3 else ('reversed' if z['n_robust_reversed']==3 else 'not_unanimously_resolved'))
        z['all_three_centre_outcome']=('preserved' if z['n_centre_preserved']==3 else ('reversed' if z['n_centre_reversed']==3 else 'not_unanimously_resolved'))
        z['exclusion_reasons']=';'.join(sorted({s for v in d.exclusion_reasons for s in str(v).split(';') if s}))
        z['flags']=';'.join(sorted({s for v in d['flags'] for s in str(v).split(';') if s}))
        ss.append(z)
    sd=pd.DataFrame(ss);sd.to_csv(HERE/'stage_pairwise_summary.csv',index=False)
    methods={
      'scope':'All distinct-gene target/EVE, target/target, target/run and target/odd event-pairs with original episode identities; counterpart pairrule-only screen stored separately.',
      'same_cohort':'Each target uses its own matched EVE temporal channel; shared bootstrap indices preserve pairing. No merged global EVE clock.',
      'independent_cohorts':'Cross-target and target/run/odd subtract independently generated bootstrap estimates on assigned ordinal stage coordinates. They remain conditional on comparable staging across cohorts.',
      'stage_availability':'Existing equally spaced stage-bin convention: observed peak before stage_index-0.5 is before current substage; within [index-0.5,index+0.5] is within current substage; beyond is later and excluded. This convention is not an estimated staging uncertainty distribution. Both CI upper bounds before current bin flagged separately.',
      'uncertainty':'20,000 biological within-stage embryo bootstrap draws, original smoothing applied to each. Pointwise percentile95% intervals conditional on feature recovery. Missing peaks are neither ties nor zero. No multiplicity claim.',
      'high_validity_flag':'At least95% bootstrap peak recovery is an explicitly pragmatic confidence flag, not a formal significance criterion; all valid fractions retained.',
      'spatial':'Original independent embryos1,2,3. Original peak extraction before posterior masking. No replacement peak found after clipping. Split EVE midpoint used only as domain centre; robust order additionally requires nonoverlapping daughter-peak spans.',
      'smoothing':'Exactly prior v04 original temporal and spatial settings, including disclosed default SVB extension. No new parameter tuning.',
      'classification':'Temporal95%contrast must exclude0 for preserved/reversed designation. Centre result and side-of-entire-lobe-span result retained separately. Unresolved and excluded records remain auditable.',
      'independence':'Pairs share genes, cycles and embryos. Counts do not represent independent tests or biological replicates.',
      'peak_identity':'Original target peak identities retained. All additional smoothed maxima listed separately and not forced into a temporal correspondence.',
    }
    (HERE/'methods.json').write_text(json.dumps(methods,indent=2))
    print(sd[(sd.stage=='Stage4_3')&(sd.all_three_robust_outcome!='not_unanimously_resolved')][['event_a','event_b','comparison_kind','all_three_robust_outcome','temporal_ci95_low','temporal_ci95_high']].to_string(index=False))

if __name__=='__main__':main()
