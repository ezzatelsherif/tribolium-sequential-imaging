"""Class-neutral peak-reference screen, original smoothing unchanged.
All pairwise eligible EVE/run/odd events are inventoried, including cross-cycle
and unresolved comparisons. Only peaks or explicit split-domain centres enter.
"""
from pathlib import Path
import itertools,json,shutil
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ORIGINAL=ROOT/'figure7_ver04_analysis'/'pairrule'
OLD=HERE/'source_inputs' if (HERE/'source_inputs').exists() else ORIGINAL
PROFILE_SOURCE=ORIGINAL if ORIGINAL.exists() else HERE
STAGE_INDEX={'Stage3_3':8,'Stage4_1':9,'Stage4_3':11,'Stage5_2':13}
MIN_VALID=.95
GENE_ORDER={'eve':0,'run':1,'odd':2}


def truth(v):return bool(v) if not pd.isna(v) else False

def canonical_tables():
    old=pd.read_csv(OLD/'spatial_features.csv');eve=old[(old.gene=='eve')&(old.event=='peak')]
    gene=pd.read_csv(OLD/'spatial_features_gene_episode.csv');gene=gene[(gene.event=='peak')&(gene.boundary_definition=='outer_episode_envelope')]
    masks=[];lookup={}
    for (stage,emb),z in eve.groupby(['stage','embryo'],sort=False):
        cycles=sorted(z.cycle.astype(int).unique(),reverse=True);selected=cycles[:3];complete=len(selected)==3
        verified_lo=float(z[z.cycle.isin(selected)].window_lo.min());requested_lo=verified_lo if complete else np.nan
        for cyc in cycles:
            r=z[z.cycle==cyc].iloc[0];rank=cycles.index(cyc)+1
            row=dict(stage=stage,embryo=int(emb),cycle=cyc,posterior_rank=rank,in_first_three=rank<=3,in_posterior_mask=rank<=3,
                ap_lo=float(r.window_lo),ap_hi=float(r.window_hi),window_lo=float(r.window_lo),window_hi=float(r.window_hi),
                mask_ap_lo=requested_lo,mask_ap_hi=1.,verified_ap_lo=verified_lo,verified_ap_hi=1.,mask_complete=complete,
                posterior_developing=truth(r.posterior_developing_cycle),identity_status='source_annotated_family',
                provenance='Original reconstructed EVE episode windows; split daughter maxima remain one domain family.',
                mask_note='' if complete else 'Only two posterior families source-identified; requested third-family boundary unresolved.')
            masks.append(row);lookup[(stage,int(emb),cyc)]=row
        if not complete:
            masks.append(dict(stage=stage,embryo=int(emb),cycle=np.nan,posterior_rank=3,in_first_three=True,in_posterior_mask=True,
                ap_lo=np.nan,ap_hi=verified_lo,window_lo=np.nan,window_hi=verified_lo,
                mask_ap_lo=np.nan,mask_ap_hi=1.,verified_ap_lo=verified_lo,verified_ap_hi=1.,mask_complete=False,
                posterior_developing=False,identity_status='unresolved_unannotated_third_family',
                provenance='Original notebook cropped the anterior region and did not annotate cycle1; full-profile audit cannot assign a consistent third family.',
                mask_note='Rank reserved; no invented cycle identity, EVE centre, or temporal correspondence.'))
    pd.DataFrame(masks).to_csv(HERE/'posterior_mask.csv',index=False)
    refs=[];sensitivity=[]
    profiles=pd.read_csv(PROFILE_SOURCE/'spatial_profiles.csv')
    for r in gene.itertuples(index=False):
        cyc=int(r.cycle);m=lookup[(r.stage,int(r.embryo),cyc)];pos=float(r.estimate) if np.isfinite(r.estimate) else np.nan
        principal=pos;lo=hi=pos;kind='peak';split=False;positions=[];flags=[]
        oldflags=str(r.flags) if not pd.isna(r.flags) else ''
        flags.extend(f for f in oldflags.split(';') if f and not any(w in f for w in ['halfheight','boundary']))
        resolved=bool(np.isfinite(pos) and truth(r.valid_estimate))
        if r.gene=='eve':
            er=eve[(eve.stage==r.stage)&(eve.embryo==r.embryo)&(eve.cycle==cyc)].iloc[0]
            pp_profile=profiles[(profiles.stage==r.stage)&(profiles.embryo==r.embryo)&(profiles.gene=='eve')]
            px=pp_profile.x_ap.to_numpy();py=pp_profile.normalized_intensity.to_numpy();pl,ph=float(er.window_lo),float(er.window_hi)
            take=(px>pl)&(px<ph);ppx=np.r_[pl,px[take],ph];ppy=np.r_[np.interp(pl,px,py),py[take],np.interp(ph,px,py)]
            selected_peaks,_=find_peaks(ppy,prominence=.10*np.ptp(ppy))
            positions=sorted(ppx[selected_peaks].tolist());split=bool(len(positions)>1 and np.ptp(positions)>=.04)
            if split:
                if len(positions)==2 and positions[1]-positions[0]>=.04:
                    lo,hi=positions;pos=(lo+hi)/2;kind='split_domain_centre';flags.append('midpoint_not_observed_peak')
                else:
                    pos=np.nan;lo=min(positions) if positions else np.nan;hi=max(positions) if positions else np.nan
                    kind='unresolved_split_domain';resolved=False;flags.append('not_exactly_two_resolved_daughters')
        tenpos,tenlo,tenhi,tenkind=pos,lo,hi,kind
        fifteenpos,fifteenlo,fifteenhi,fifteenkind=pos,lo,hi,kind
        criterion_sensitive=False
        if r.gene=='eve':
            profile=profiles[(profiles.stage==r.stage)&(profiles.embryo==r.embryo)&(profiles.gene=='eve')]
            x=profile.x_ap.to_numpy();y=profile.normalized_intensity.to_numpy();a,b=float(er.window_lo),float(er.window_hi)
            inside=(x>a)&(x<b);xx=np.r_[a,x[inside],b];yy=np.r_[np.interp(a,x,y),y[inside],np.interp(b,x,y)]
            for fraction in [.10,.15]:
                peaks,_=find_peaks(yy,prominence=fraction*np.ptp(yy));pp=xx[peaks]
                split_here=len(pp)>1 and np.ptp(pp)>=.04
                if split_here and len(pp)==2:
                    cp=float(np.mean(pp));cl=float(pp.min());ch=float(pp.max());ck='split_domain_centre'
                elif split_here:
                    cp=np.nan;cl=float(pp.min());ch=float(pp.max());ck='unresolved_split_domain'
                else:cp=principal;cl=principal;ch=principal;ck='peak'
                sensitivity.append(dict(stage=r.stage,embryo=int(r.embryo),gene='eve',cycle=cyc,event_id=f'eve_cycle{cyc}',
                    prominence_fraction=fraction,min_split_separation_ap=.04,n_prominent_peaks=len(pp),daughter_peaks=json.dumps(pp.tolist()),
                    position=cp,lobe_lo=cl,lobe_hi=ch,reference_kind=ck,in_posterior_mask=bool(m['in_first_three'])))
                if fraction==.10:tenpos,tenlo,tenhi,tenkind=cp,cl,ch,ck
                else:fifteenpos,fifteenlo,fifteenhi,fifteenkind=cp,cl,ch,ck
            criterion_sensitive=bool(fifteenkind!=kind or not np.isclose(fifteenpos,pos,equal_nan=True))
            if criterion_sensitive:flags.append('split_detection_sensitive_10_vs_15pct_prominence')
        inmask=bool(m['in_first_three'])
        inside_verified=bool(np.isfinite(pos) and pos>=m['verified_ap_lo'] and pos<=1.)
        if inmask and np.isfinite(pos) and not inside_verified:flags.append('reference_outside_verified_spatial_mask')
        if m['posterior_developing']:flags.append('posterior_developing_family')
        refs.append(dict(stage=r.stage,embryo=int(r.embryo),gene=r.gene,cycle=cyc,event_id=f'{r.gene}_cycle{cyc}',
            posterior_rank=m['posterior_rank'],position=pos,lobe_lo=lo,lobe_hi=hi,span_lo=lo,span_hi=hi,
            reference_kind=kind,principal_peak=principal,split_flag=split,daughter_peaks=json.dumps(positions),
            split_criterion_sensitive=criterion_sensitive,split_prominence_fraction=.10,tenpct_position=tenpos,tenpct_lobe_lo=tenlo,tenpct_lobe_hi=tenhi,tenpct_reference_kind=tenkind,
            fifteenpct_position=fifteenpos,fifteenpct_lobe_lo=fifteenlo,fifteenpct_lobe_hi=fifteenhi,fifteenpct_reference_kind=fifteenkind,
            spatial_resolved=resolved,in_posterior_mask=inmask,in_verified_spatial_mask=inside_verified,
            posterior_developing=m['posterior_developing'],mask_complete=m['mask_complete'],
            mask_ap_lo=m['mask_ap_lo'],mask_ap_hi=1.,verified_ap_lo=m['verified_ap_lo'],
            family_window_lo=m['window_lo'],family_window_hi=m['window_hi'],
            gene_window_lo=r.gene_window_lo,gene_window_hi=r.gene_window_hi,
            flags=';'.join(sorted(set(flags)))))
    pd.DataFrame(refs).to_csv(HERE/'spatial_peak_references.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(HERE/'spatial_eve_split_sensitivity.csv',index=False)
    t=pd.read_csv(OLD/'temporal_features.csv');t=t[t.event=='peak'].copy();t['event_id']=t.gene+'_cycle'+t.cycle.astype(int).astype(str)
    t['reference_kind']='peak';t['position']=t.estimate;t['smoothed_fit_conditional']=True
    t.to_csv(HERE/'temporal_peak_references.csv',index=False)
    z=np.load(OLD/'temporal_feature_bootstraps.npz');draws={r.event_id:z[r.event_id+'_peak'] for r in t.itertuples()}
    np.savez_compressed(HERE/'temporal_peak_bootstraps.npz',**draws)
    for name in ['temporal_profiles.csv','spatial_profiles.csv','temporal_raw_observations.csv','spatial_raw_observations.csv']:
        if (PROFILE_SOURCE/name).resolve()!=(HERE/name).resolve():shutil.copyfile(PROFILE_SOURCE/name,HERE/name)
    return pd.DataFrame(refs),t,draws


def statistics(d):
    valid=np.isfinite(d[1:]);v=d[1:][valid];ci=np.quantile(v,[.025,.975]) if len(v) and np.isfinite(d[0]) else [np.nan,np.nan]
    return dict(estimate=float(d[0]),ci95_low=float(ci[0]),ci95_high=float(ci[1]),valid_bootstrap_n=int(valid.sum()),
        valid_bootstrap_fraction=float(valid.mean()),bootstrap_fraction_positive=float(np.mean(v>0)) if len(v) else np.nan,
        bootstrap_fraction_equal=float(np.mean(v==0)) if len(v) else np.nan)


def pair_screen(sp,t,draws):
    temporal_lookup=t.set_index('event_id');temporal_pairs={}
    ids=sorted(draws,key=lambda k:(int(k.split('cycle')[1]),GENE_ORDER[k.split('_')[0]]))
    for a,b in itertools.combinations(ids,2):
        temporal_pairs[(a,b)]=statistics(draws[b]-draws[a])
    rows=[]
    for (stage,emb),z in sp[sp.in_posterior_mask].groupby(['stage','embryo'],sort=False):
        z=z.sort_values(['cycle','gene'],key=lambda s:s.map(GENE_ORDER) if s.name=='gene' else s)
        for a,b in itertools.combinations(list(z.itertuples(index=False)),2):
            meta=temporal_pairs[(a.event_id,b.event_id)];delta=b.position-a.position
            t_sign=1 if meta['ci95_low']>0 else (-1 if meta['ci95_high']<0 else 0)
            s_sign=int(np.sign(delta)) if np.isfinite(delta) else 0
            overlap=not(b.lobe_lo>a.lobe_hi or a.lobe_lo>b.lobe_hi) if np.isfinite([a.lobe_lo,a.lobe_hi,b.lobe_lo,b.lobe_hi]).all() else True
            span_sign=1 if np.isfinite([b.lobe_lo,a.lobe_hi]).all() and b.lobe_lo>a.lobe_hi else (-1 if np.isfinite([a.lobe_lo,b.lobe_hi]).all() and a.lobe_lo>b.lobe_hi else 0)
            conditional='preserved' if t_sign and s_sign==t_sign else ('reversed' if t_sign and s_sign==-t_sign else 'unresolved')
            stage_end=STAGE_INDEX[stage];bin_lo=stage_end-.5;bin_hi=stage_end+.5;ra=temporal_lookup.loc[a.event_id];rb=temporal_lookup.loc[b.event_id]
            def availability(value):
                if not np.isfinite(value):return 'unresolved'
                if value<bin_lo:return 'before_current_substage'
                if value<=bin_hi:return 'within_current_substage'
                return 'later_substage'
            av_a,av_b=availability(ra.estimate),availability(rb.estimate)
            valid_draws=np.isfinite(draws[a.event_id][1:])&np.isfinite(draws[b.event_id][1:])
            before=valid_draws&(draws[a.event_id][1:]<bin_lo)&(draws[b.event_id][1:]<bin_lo)
            by_bin_end=valid_draws&(draws[a.event_id][1:]<=bin_hi)&(draws[b.event_id][1:]<=bin_hi)
            ci_before=bool(np.isfinite(ra.ci95_high) and np.isfinite(rb.ci95_high) and max(ra.ci95_high,rb.ci95_high)<bin_lo)
            ten_span_sign=1 if np.isfinite([b.tenpct_lobe_lo,a.tenpct_lobe_hi]).all() and b.tenpct_lobe_lo>a.tenpct_lobe_hi else (-1 if np.isfinite([a.tenpct_lobe_lo,b.tenpct_lobe_hi]).all() and a.tenpct_lobe_lo>b.tenpct_lobe_hi else 0)
            fifteen_span_sign=1 if np.isfinite([b.fifteenpct_lobe_lo,a.fifteenpct_lobe_hi]).all() and b.fifteenpct_lobe_lo>a.fifteenpct_lobe_hi else (-1 if np.isfinite([a.fifteenpct_lobe_lo,b.fifteenpct_lobe_hi]).all() and a.fifteenpct_lobe_lo>b.fifteenpct_lobe_hi else 0)
            split_stable=bool(span_sign==fifteen_span_sign and span_sign!=0)
            completed=np.isfinite(draws[a.event_id][1:])&np.isfinite(draws[b.event_id][1:])&(draws[a.event_id][1:]<=stage_end)&(draws[b.event_id][1:]<=stage_end)
            observed_completed=bool(np.isfinite(ra.estimate) and np.isfinite(rb.estimate) and max(ra.estimate,rb.estimate)<=stage_end)
            ci_completed=bool(np.isfinite(ra.ci95_high) and np.isfinite(rb.ci95_high) and max(ra.ci95_high,rb.ci95_high)<=stage_end)
            exclusions=[];flags=[]
            if not(a.spatial_resolved and b.spatial_resolved):exclusions.append('spatial_reference_unresolved')
            if (np.isfinite(a.position) and not a.in_verified_spatial_mask) or (np.isfinite(b.position) and not b.in_verified_spatial_mask):exclusions.append('outside_verified_spatial_mask')
            if a.posterior_developing or b.posterior_developing:exclusions.append('developing_family')
            if 'later_substage' in [av_a,av_b]:exclusions.append('temporal_peak_in_later_substage')
            if 'unresolved' in [av_a,av_b]:exclusions.append('temporal_peak_unresolved')
            if 'within_current_substage' in [av_a,av_b]:flags.append('temporal_peak_within_current_substage')
            if not ci_before:flags.append('temporal_history_before_current_substage_not_secured_by_ci')
            high_validity=meta['valid_bootstrap_fraction']>=MIN_VALID
            if not high_validity:flags.append('conditional_interval_with_lower_peak_validity')
            if overlap and(a.split_flag or b.split_flag):flags.append('order_depends_on_split_domain_centre')
            if span_sign!=fifteen_span_sign:flags.append('comparison_sensitive_to_eve_split_prominence_criterion')
            if not a.mask_complete:flags.append('third_posterior_family_unresolved_but_verified_subset_used')
            eligible=not exclusions
            if not t_sign:
                conditional='unresolved_temporal_order';robust='unresolved_temporal_order';robust15='unresolved_temporal_order'
            else:
                conditional='preserved' if s_sign==t_sign else ('reversed' if s_sign==-t_sign else 'unresolved_spatial_order')
                robust='preserved' if span_sign==t_sign else ('reversed' if span_sign==-t_sign else 'within_split_domain_or_tied')
                robust15='preserved' if fifteen_span_sign==t_sign else ('reversed' if fifteen_span_sign==-t_sign else 'within_split_domain_or_tied')
            if not eligible:conditional=robust=robust15='excluded'
            strong=bool(eligible and high_validity and ci_before and robust in ['preserved','reversed'] and robust15==robust)
            reasons=exclusions.copy()
            if robust=='unresolved_temporal_order':reasons.append('temporal_order_ci_includes_zero_or_unresolved')
            if robust=='within_split_domain_or_tied':reasons.append('within_split_domain_or_tied')
            crosscycle=a.cycle!=b.cycle
            windows_separated=bool(ra.window_hi_index<=rb.window_lo_index or rb.window_hi_index<=ra.window_lo_index)
            rows.append(dict(stage=stage,embryo=emb,event_a=a.event_id,event_b=b.event_id,gene_a=a.gene,gene_b=b.gene,
                cycle_a=a.cycle,cycle_b=b.cycle,posterior_rank_a=a.posterior_rank,posterior_rank_b=b.posterior_rank,
                reference_kind_a=a.reference_kind,reference_kind_b=b.reference_kind,position_a=a.position,position_b=b.position,
                lobe_a_lo=a.lobe_lo,lobe_a_hi=a.lobe_hi,lobe_b_lo=b.lobe_lo,lobe_b_hi=b.lobe_hi,
                temporal_delta_b_minus_a=meta['estimate'],temporal_ci95_low=meta['ci95_low'],temporal_ci95_high=meta['ci95_high'],
                temporal_valid_bootstrap_n=meta['valid_bootstrap_n'],temporal_valid_bootstrap_fraction=meta['valid_bootstrap_fraction'],
                temporal_fraction_positive=meta['bootstrap_fraction_positive'],spatial_delta_b_minus_a=delta,
                spatial_span_gap_b_after_a=b.lobe_lo-a.lobe_hi,spatial_span_gap_a_after_b=a.lobe_lo-b.lobe_hi,
                temporal_order_sign=t_sign,spatial_centre_order_sign=s_sign,spatial_span_order_sign=span_sign,
                split_lobe_span_robust=not overlap,centre_based_outcome=conditional,robust_outcome=robust,fifteenpct_robust_outcome=robust15,
                eligible=eligible,exclusion_reasons=';'.join(exclusions),qualification_flags=';'.join(flags),high_feature_validity=high_validity,
                matched_temporal_channels=True,strong_claim_eligible=strong,strong_within_cycle_eligible=strong and a.cycle==b.cycle,
                observed_temporal_peaks_completed_by_stage=observed_completed,temporal_ci_completed_by_stage=ci_completed,
                temporal_joint_completed_fraction=float(completed.mean()),temporal_stage_limit=stage_end,
                temporal_substage_bin_lo=bin_lo,temporal_substage_bin_hi=bin_hi,temporal_availability_a=av_a,temporal_availability_b=av_b,
                temporal_joint_before_current_substage_fraction=float(before.mean()),temporal_joint_by_current_substage_end_fraction=float(by_bin_end.mean()),
                temporal_ci_before_current_substage=ci_before,spatial_span_order_sign_10pct=ten_span_sign,spatial_span_order_sign_15pct=fifteen_span_sign,split_span_order_stable_10_15pct=split_stable,
                cross_cycle=crosscycle,window_separated=windows_separated,same_gene=a.gene==b.gene,
                support_role='identity_conditional_cross_cycle_inventory' if crosscycle else 'within_cycle_comparison',
                confidence_intervals_multiplicity_adjusted=False,mask_complete=a.mask_complete,
                unresolved_reasons=';'.join(reasons),flags_a=a.flags,flags_b=b.flags))
    out=pd.DataFrame(rows);out.to_csv(HERE/'pairwise_peak_screen.csv',index=False)
    summary=[]
    for (stage,a,b),d in out.groupby(['stage','event_a','event_b'],sort=False):
        summary.append(dict(stage=stage,event_a=a,event_b=b,n_embryos=len(d),
            n_centre_preserved=int((d.centre_based_outcome=='preserved').sum()),n_centre_reversed=int((d.centre_based_outcome=='reversed').sum()),
            n_robust_preserved=int((d.robust_outcome=='preserved').sum()),n_robust_reversed=int((d.robust_outcome=='reversed').sum()),
            n_unresolved=int(d.robust_outcome.isin(['unresolved_temporal_order','within_split_domain_or_tied']).sum()),
            n_excluded=int((d.robust_outcome=='excluded').sum()),n_eligible_embryos=int(d.eligible.sum()),
            n_strong_claim_eligible=int(d.strong_claim_eligible.sum()),n_strong_preserved=int(((d.robust_outcome=='preserved')&d.strong_claim_eligible).sum()),
            n_strong_reversed=int(((d.robust_outcome=='reversed')&d.strong_claim_eligible).sum()),
            all_three_robust_preserved=bool(len(d)==3 and (d.robust_outcome=='preserved').all()),
            all_three_robust_reversed=bool(len(d)==3 and (d.robust_outcome=='reversed').all()),
            all_three_strong_preserved=bool(len(d)==3 and ((d.robust_outcome=='preserved')&d.strong_claim_eligible).all()),
            all_three_strong_reversed=bool(len(d)==3 and ((d.robust_outcome=='reversed')&d.strong_claim_eligible).all()),
            all_three_centre_preserved=bool(len(d)==3 and (d.centre_based_outcome=='preserved').all()),
            temporal_delta_b_minus_a=d.iloc[0].temporal_delta_b_minus_a,temporal_ci95_low=d.iloc[0].temporal_ci95_low,
            temporal_ci95_high=d.iloc[0].temporal_ci95_high,temporal_valid_fraction=d.iloc[0].temporal_valid_bootstrap_fraction,
            cross_cycle=bool(d.iloc[0].cross_cycle),window_separated=bool(d.iloc[0].window_separated),
            same_gene=bool(d.iloc[0].same_gene),support_role=d.iloc[0].support_role,
            unresolved_reasons=';'.join(sorted(set(';'.join(d.unresolved_reasons.fillna('')).split(';'))- {''})),
            qualification_flags=';'.join(sorted(set(';'.join(d.qualification_flags.fillna('')).split(';'))- {''}))))
    pd.DataFrame(summary).to_csv(HERE/'pairwise_stage_summary.csv',index=False)
    tp=[]
    for (a,b),v in temporal_pairs.items():tp.append(dict(event_a=a,event_b=b,**v))
    pd.DataFrame(tp).to_csv(HERE/'all_temporal_pairwise_contrasts.csv',index=False)
    return out,pd.DataFrame(summary)

if __name__=='__main__':
    sp,t,draws=canonical_tables();out,s=pair_screen(sp,t,draws)
    print('Spatial references:',len(sp),'pair rows:',len(out),'stage-pair summaries:',len(s))
    print(s[(~s.cross_cycle)&s.all_three_robust_preserved][['stage','event_a','event_b','temporal_valid_fraction']].to_string(index=False))
