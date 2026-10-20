"""GT/MLPT direct features after the original handover smoothing.
No WDL calculation. Curated episode windows remain fixed across raw/smoothed
comparisons; a selected maximum must be an interior local maximum to be a peak.
Bootstrap stage means are re-smoothed individually before feature extraction.
"""
from pathlib import Path
import json, sys
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parents[1]/'all_gene_direct_analysis'
sys.path.insert(0,str(ROOT.parent))
from smoothing import temporal,spatial,spatial_parameters,window_length
LABELS=[f'{c}.{s}' for c in range(1,9) for s in range(1,4)]
EPISODES={'GT':(6,15),'MLPT':(2,12)}
PEAKS={'GT':[('pulse1',6,9),('pulse2',9,11)],'MLPT':[('first',2,12)]}

def stats(v):
 v=np.asarray(v);good=np.isfinite(v[1:]);b=v[1:][good];ci=np.quantile(b,[.025,.975]) if good.any() and np.isfinite(v[0]) else [np.nan,np.nan]
 return dict(estimate=v[0],ci95_low=ci[0],ci95_high=ci[1],valid_bootstrap_n=int(good.sum()),valid_bootstrap_fraction=float(good.mean()),valid_fraction=float(good.mean()))

def measurements(x,y,g):
 n=len(y);r=np.arange(n);lo,hi=EPISODES[g];use=(x>=lo)&(x<=hi);base=y.min(1);maximum=y[:,use].max(1);level=base+.5*(maximum-base)
 above=y>=level[:,None];ccrise=(~above[:,:-1])&above[:,1:];ccfall=above[:,:-1]&(~above[:,1:]);region=(x[:-1]>=lo)&(x[1:]<=hi)
 out={'baseline':base,'episode_height':maximum,'threshold':level}
 for ev,cc in [('rise',ccrise),('fall',ccfall)]:
  cc&=region[None,:];ii=cc.argmax(1) if ev=='rise' else cc.shape[1]-1-cc[:,::-1].argmax(1)
  with np.errstate(divide='ignore',invalid='ignore'):v=x[ii]+(level-y[r,ii])*(x[ii+1]-x[ii])/(y[r,ii+1]-y[r,ii])
  v[~cc.any(1)]=np.nan;out[ev]=v
 for ev,a,b in PEAKS[g]:
  ids=np.flatnonzero((x>=a)&(x<=b));ii=ids[np.argmax(y[:,ids],1)];ok=(ii>0)&(ii<len(x)-1)
  ok&=(y[r,ii]>=y[r,np.maximum(0,ii-1)])&(y[r,ii]>=y[r,np.minimum(len(x)-1,ii+1)])
  v=x[ii].astype(float);v[~ok]=np.nan;h=y[r,ii].copy();h[~ok]=np.nan
  out[ev+'_peak']=v;out[ev+'_height']=h
 if g=='GT':
  a=out['pulse1_peak'];b=out['pulse2_peak'];mask=(x[None,:]>=a[:,None])&(x[None,:]<=b[:,None]);yy=np.where(mask,y,np.inf);i=yy.argmin(1);v=yy[r,i];v[~np.isfinite(v)]=np.nan
  out['interpeak_trough']=v;out['interpeak_trough_x']=np.where(np.isfinite(v),x[i],np.nan)
  out['gt2_over_gt1_above_series_min']=(out['pulse2_height']-base)/(out['pulse1_height']-base)
  out['valley_fraction_of_smaller_peak']=(v-base)/(np.minimum(out['pulse1_height'],out['pulse2_height'])-base)
 return out

def cut(x,y,lo,hi):
 take=(x>lo)&(x<hi);return np.r_[lo,x[take],hi],np.r_[np.interp(lo,x,y),y[take],np.interp(hi,x,y)]
def cross(x,y,q,side):
 ix=np.flatnonzero(y>=q)
 if not len(ix):return np.nan,True
 k=ix[0] if side=='rise' else ix[-1]
 if (side=='rise' and k==0) or (side=='fall' and k==len(x)-1):return np.nan,True
 a,b=(k-1,k) if side=='rise' else (k,k+1)
 return x[a]+(q-y[a])*(x[b]-x[a])/(y[b]-y[a]),False

def run_temporal():
 curves=np.load(BASE/'temporal_bootstrap_curves.npz');all_features={};raw_features={};tr=[];pk=[];profiles=[];raw_profiles=[];metricrows=[];boot={};rawboot={}
 for gene in ['GT','MLPT']:
  raw=curves[gene+'_'+gene];raw_features[gene]=measurements(np.arange(24),raw,gene)
  # 20,000 biological draws already resample independent embryos within each stage.
  dense=[]
  for start in range(0,len(raw),500):
   x,y=temporal(raw[start:start+500]);dense.append(y)
  y=np.concatenate(dense,axis=0);f=measurements(x,y,gene);all_features[gene]=f
  base=f['baseline'];amp=f['episode_height']-base;scaled=(y-base[:,None])/amp[:,None]
  low,high=np.quantile(scaled[1:],[.025,.975],axis=0)
  profiles.append(pd.DataFrame({'gene':gene,'substage_index':x,'normalized_mean':(y[0]-base[0])/(y[0].max()-base[0]),'episode_scaled_mean':scaled[0],'smoothed_original_scale':y[0],'ci95_low':low,'ci95_high':high}))
  # Raw dots transformed with the same original full-series normalization and
  # observed smooth episode scale, so they remain comparable with the curve.
  rr=(raw[0]-raw[0].min())/np.ptp(raw[0]);raw_profiles.append(pd.DataFrame({'gene':gene,'stage':LABELS,'substage_index':range(24),'raw_mean':raw[0],'original_normalized_mean':rr,'episode_scaled_mean':(rr-base[0])/amp[0]}))
  for event in ['rise','fall']:
   tr.append(dict(gene=gene,event=event,normalization='episode',threshold=.5,anchor_method='none',window_start=LABELS[EPISODES[gene][0]],window_end=LABELS[EPISODES[gene][1]],unit='ordinal_substage',**stats(f[event])))
  for event,lo,hi in PEAKS[gene]:
   v=f[event+'_peak'];pk.append(dict(gene=gene,event=event,anchor_method='none',window_start=LABELS[lo],window_end=LABELS[hi],unit='ordinal_substage',peak_height=f[event+'_height'][0],normalized_peak_height=(f[event+'_height'][0]-base[0])/amp[0],**stats(v)))
  for key,v in f.items():
   boot[gene+'_'+key]=v;rawboot[gene+'_'+key]=raw_features[gene][key]
   if key not in ['baseline','threshold','episode_height']:
    metricrows.append(dict(gene=gene,metric=key,processing='smoothed',**stats(v)))
    metricrows.append(dict(gene=gene,metric=key,processing='raw',**stats(raw_features[gene][key])))
 pd.DataFrame(tr).to_csv(ROOT/'temporal_gap_events.csv',index=False);pd.DataFrame(pk).to_csv(ROOT/'temporal_gap_peaks.csv',index=False)
 pd.concat(profiles).to_csv(ROOT/'temporal_profiles.csv',index=False);pd.concat(raw_profiles).to_csv(ROOT/'temporal_raw_profiles.csv',index=False)
 pd.DataFrame(metricrows).to_csv(ROOT/'temporal_raw_vs_smoothed.csv',index=False)
 contrasts=[]
 for mode,fs in [('smoothed',all_features),('raw',raw_features)]:
  for name,a,b in [('GT1_minus_MLPT_peak','pulse1_peak','first_peak'),('GT2_minus_MLPT_peak','pulse2_peak','first_peak'),('GT_minus_MLPT_rise','rise','rise'),('GT_minus_MLPT_fall','fall','fall')]:
   v=fs['GT'][a]-fs['MLPT'][b];contrasts.append(dict(processing=mode,contrast=name,unit='ordinal_substage',bootstrap_fraction_positive=float(np.mean(v[1:]>0)),**stats(v)))
 pd.DataFrame(contrasts).to_csv(ROOT/'temporal_contrasts.csv',index=False)
 np.savez_compressed(ROOT/'temporal_bootstrap_features.npz',**boot)
 return all_features

def run_spatial():
 raw=pd.read_csv(BASE/'spatial_raw_profiles.csv');cfgs=json.loads((BASE/'reconstructed_points.json').read_text())['configs']
 profiles=[];domains=[];episodes=[];allmax=[];sens=[];summary=[];comp=[]
 for cfg in cfgs:
  stage=cfg['stage']
  for embryo in [1,2,3]:
   data={}
   for gene in ['gt','mlpt']:
    z=raw[(raw.stage==stage)&(raw.embryo==embryo)&(raw.gene==gene)];x=z.x_ap.to_numpy();rr=z.raw_intensity.to_numpy();_,s=spatial(x,rr,stage,embryo,gene);rn=z.normalized_intensity.to_numpy()
    lo=.45;ids=np.flatnonzero(x>=lo);base=s.min();peak=s[ids].max();level=base+.5*(peak-base);scaled=(s-base)/(peak-base)
    profiles.append(pd.DataFrame(dict(stage=stage,embryo=embryo,gene=gene,x_ap=x,raw_intensity=rr,raw_normalized_intensity=rn,normalized_intensity=scaled,smoothed_original_scale=s,raw_episode_scaled=(rn-base)/(peak-base))))
    data[gene]=(x,s,rn,rr)
    for processing,y in [('smoothed',s),('raw',rn)]:
     b=y.min();p=y[ids].max();q=b+.5*(p-b);xx,yy=cut(x,y,lo,1);rise,rc=cross(xx,yy,q,'rise');fall,fc=cross(xx,yy,q,'fall')
     rec=dict(stage=stage,embryo=embryo,gene=gene,threshold_fraction=.5,window_lo=lo,window_hi=1.,peak=x[ids[np.argmax(y[ids])]],rise=rise,fall=fall,rise_censored=rc,fall_censored=fc,excluded_neural=False,threshold_level=q,baseline=b,episode_height=p,processing=processing)
     comp.append(rec)
     if processing=='smoothed':domains.append(rec)
    points=cfg['points'][str(embryo)][gene]
    for i in range(1,len(points),2):
     n=i//2+1;lo=points[i-1];hi=points[i+1] if i+1<len(points) else 1.;ids=np.flatnonzero((x>=lo)&(x<=hi))
     for processing,y in [('smoothed',s),('raw',rn)]:
      k=int(ids[np.argmax(y[ids])]);resolved=0<k<len(x)-1 and y[k]>=y[k-1] and y[k]>=y[k+1];event='gt_peak'+str(n) if gene=='gt' else 'mlpt_main'
      rec=dict(stage=stage,embryo=embryo,gene=gene,group='gap',event=event,cycle=np.nan,complete_coarse_window=i+1<len(points),excluded_neural=False,peak=x[k] if resolved else np.nan,candidate_maximum_x=x[k],peak_height=y[k] if resolved else np.nan,window_lo=lo,window_hi=hi,peak_at_window_edge=k in [ids[0],ids[-1]],is_local_maximum=resolved,processing=processing,normalized_peak_height=(y[k]-y.min())/(y.max()-y.min()) if resolved else np.nan)
      if processing=='smoothed':episodes.append(rec)
      summary.append(rec)
   x,y,rn,rr=data['gt'];postlo=cfg['points'][str(embryo)]['gt'][0];ids=np.flatnonzero(x>=postlo);xx=x[ids];yy=y[ids];pp,pr=find_peaks(yy,prominence=0);w=peak_widths(yy,pp,rel_height=.5);wl=np.interp(w[2],np.arange(len(xx)),xx);wr=np.interp(w[3],np.arange(len(xx)),xx);span=np.ptp(yy)
   for k,prom,left,right in zip(pp,pr['prominences'],wl,wr):allmax.append(dict(stage=stage,embryo=embryo,gene='gt',peak=x[ids[k]],height=yy[k],prominence_posterior_range=prom/span,half_prominence_width_ap=right-left))
   for frac in [.05,.10,.15,.20]:
    for width in [.02,.04]:
     good=np.flatnonzero((pr['prominences']>=frac*span)&((wr-wl)>=width));kept=[]
     for j in good[np.argsort(pr['prominences'][good])[::-1]]:
      if all(abs(xx[pp[j]]-xx[pp[v]])>=.04 for v in kept):kept.append(j)
     kept=sorted(kept,key=lambda j:xx[pp[j]])
     sens.append(dict(stage=stage,embryo=embryo,min_prominence_posterior_range=frac,min_width_ap=width,n_broad_peaks=len(kept),broad_peak_positions=json.dumps([float(xx[pp[j]]) for j in kept])))
 pd.concat(profiles).to_csv(ROOT/'spatial_profiles.csv',index=False);pd.DataFrame(domains).to_csv(ROOT/'spatial_domain_metrics.csv',index=False);pd.DataFrame(episodes).to_csv(ROOT/'spatial_episode_metrics.csv',index=False)
 pd.DataFrame(comp).to_csv(ROOT/'spatial_halfheight_raw_vs_smoothed.csv',index=False);pd.DataFrame(summary).to_csv(ROOT/'spatial_peaks_raw_vs_smoothed.csv',index=False);pd.DataFrame(allmax).to_csv(ROOT/'spatial_all_gt_local_maxima.csv',index=False);pd.DataFrame(sens).to_csv(ROOT/'spatial_broad_peak_sensitivity.csv',index=False)
 # Each stage is a different set of embryos; no tracking of individual embryos over time.
 peaks=pd.DataFrame(summary);prof=pd.concat(profiles);rows=[]
 for (stage,embryo,processing),d in peaks.groupby(['stage','embryo','processing']):
  e=d.set_index('event');p=prof[(prof.stage==stage)&(prof.embryo==embryo)&(prof.gene=='gt')];y=p.smoothed_original_scale if processing=='smoothed' else p.raw_normalized_intensity;x=p.x_ap.to_numpy();y=y.to_numpy();r=p.raw_intensity.to_numpy();out=dict(stage=stage,embryo=embryo,processing=processing,gt1_x=e.loc['gt_peak1','peak'],mlpt_peak_x=e.loc['mlpt_main','peak'],gt2_x=np.nan,gt2_over_gt1_above_full_min=np.nan,gt2_over_gt1_reconstructed_intensity=np.nan,mlpt_between_gt_peaks=None)
  if 'gt_peak2' in e.index and np.isfinite(e.loc['gt_peak2','peak']):
   a,b=e.loc['gt_peak1','peak'],e.loc['gt_peak2','peak'];ha=np.interp(a,x,y);hb=np.interp(b,x,y);u=(x>=a)&(x<=b);v=y[u].min();out.update(gt2_x=b,gt2_over_gt1_above_full_min=(hb-y.min())/(ha-y.min()),gt2_over_gt1_reconstructed_intensity=(hb*np.ptp(r)+r.min())/(ha*np.ptp(r)+r.min()),mlpt_between_gt_peaks=bool(a<out['mlpt_peak_x']<b),valley_fraction_of_smaller_peak=(v-y.min())/(min(ha,hb)-y.min()),dominant_lobe=1 if ha>=hb else 2)
  rows.append(out)
 pd.DataFrame(rows).to_csv(ROOT/'spatial_two_lobe_comparison.csv',index=False)
 m=peaks[peaks.processing=='raw'].merge(peaks[peaks.processing=='smoothed'],on=['stage','embryo','gene','event'],suffixes=('_raw','_smoothed'));m['shift_percent_ap']=100*(m.peak_smoothed-m.peak_raw)
 m[['stage','embryo','gene','event','peak_raw','peak_smoothed','shift_percent_ap']].to_csv(ROOT/'spatial_peak_displacement.csv',index=False)
 rows=[]
 for (s,e),a in pd.DataFrame(domains).groupby(['stage','embryo']):
  a=a.set_index('gene')
  for ev in ['rise','fall']:rows.append({'stage':s,'embryo':e,'event':ev,'GT_minus_MLPT_percent_AP':100*(a.loc['gt',ev]-a.loc['mlpt',ev])})
 pd.DataFrame(rows).to_csv(ROOT/'spatial_halfheight_contrasts.csv',index=False)
 extra=[]
 for r in pd.DataFrame(allmax).query('prominence_posterior_range >= .10 and half_prominence_width_ap >= .04').itertuples():
  assigned=pd.DataFrame(episodes).query('stage == @r.stage and embryo == @r.embryo and gene == \"gt\"').peak
  if all(abs(v-r.peak)>.04 for v in assigned):extra.append({k:v for k,v in r._asdict().items() if k!='Index'})
 pd.DataFrame(extra).to_csv(ROOT/'spatial_additional_unassigned_gt_peaks.csv',index=False)

if __name__=='__main__':
 run_temporal();run_spatial()
 methods={'temporal_smoothing':'Original full-series min-max normalization, PCHIP from24 ordered stage means to800 points; SG window81 polynomial3 mode interp. Each biological bootstrap curve independently smoothed.', 'spatial_smoothing':'Original full-profile minmax normalization, exact per-stage/embryo/gene fractional SG windows before cropping; no spatial resampling.', 'halfheight':'Full smoothed series/profile minimum +0.5*(selected smoothed episode maximum−minimum). No clipping of SG overshoot. Raw comparator applies corresponding rule to raw.', 'temporal_windows':EPISODES,'temporal_peak_windows':PEAKS,'spatial_windows':'Same curated v03 GT/MLPT lobe windows, full outer posterior domain x=.45..1. A maximum must be a true interior local maximum, otherwise unresolved. Additional smoothed maxima outside assigned windows audited separately.', 'bootstrap':'20,000 within-stage biological bootstrap curves from all_gene_direct_analysis/temporal_bootstrap_curves.npz. GT and MLPT are separate cohorts. CIs conditional on stage assignment, smoothing method and valid feature.', 'stage_caution':'Substage coordinates are ordinal, not calibrated times. Dense PCHIP does not create independent time observations; smooth peak resolution may exceed underlying sampling resolution.', 'spatial_replicates':'Original three independent embryos each stage. Different embryos across stages; pixels not independent replicates.'}
 (ROOT/'methods.json').write_text(json.dumps(methods,indent=2))
 print(pd.read_csv(ROOT/'temporal_contrasts.csv').to_string(index=False));print(pd.read_csv(ROOT/'spatial_two_lobe_comparison.csv').query("processing=='smoothed'").to_string(index=False))
