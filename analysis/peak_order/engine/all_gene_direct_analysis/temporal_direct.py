"""Unsmoothed temporal events with paired eve registration and embryo resampling.
Run: python temporal_direct.py --bootstraps 20000
This is fixed-embryo pseudotime, not longitudinal trajectories or calibrated minutes.
"""
from pathlib import Path
import argparse,json,hashlib,warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'data'/'Temporal'
GAPS=['KR','SVB','MLPT','GT']; COLORS={'KR':'#1565a3','SVB':'#7042a5','MLPT':'#248342','GT':'#c27d0c','eve':'#923640','run':'#225ca2','odd':'#469247'}
LABELS=[f'{c}.{s}' for c in range(1,9) for s in range(1,4)]
WINDOWS={'KR':('1.1','4.3'),'SVB':('1.1','4.3'),'MLPT':('1.3','5.1'),'GT':('3.1','6.1')}

def summary(v):
 v=np.asarray(v);good=np.isfinite(v[1:]);b=v[1:][good];obs=float(v[0]) if np.isfinite(v[0]) else None
 # Conditional intervals are reported only when observed statistic exists; valid fraction is always exposed.
 ci=np.quantile(b,[.025,.975]) if len(b) and obs is not None else [np.nan,np.nan]
 return dict(estimate=obs,ci95_low=float(ci[0]),ci95_high=float(ci[1]),valid_bootstrap_n=int(len(b)),valid_bootstrap_fraction=float(good.mean()))

def brackets(t):
 if not np.isfinite(t):return ('','')
 lo=int(np.floor(t));hi=int(np.ceil(t));return (LABELS[max(0,min(23,lo))],LABELS[max(0,min(23,hi))])

def load_sources():
 cohorts={};audit={};tidy=[]
 for g in ['run','odd']+GAPS+['HB']:
  p=SOURCE/f'{g}.xlsx';d=pd.read_excel(p,sheet_name='Sheet1');base=pd.read_excel(p,sheet_name='BASE');base.columns=base.columns.str.strip()
  assert [f'{s:.1f}' for s in d.Stage]==LABELS
  dat={ch:np.stack([(d[f'{ch}{k}A']-d[f'{ch}{k}B']).to_numpy(float) for k in range(1,6)],1) for ch in ['eve','mrna']}
  aa={}
  for ch,a in dat.items():
   cache=base[[f'{ch}{k}' for k in range(1,6)]].to_numpy(float);mask=np.isfinite(a)
   aa[ch]={'observed':int(mask.sum()),'missing':int((~mask).sum()),'BASE_matches_valid_A_minus_B':bool(np.allclose(a[mask],cache[mask])),'BASE_missing_slots_zero':bool(np.all(cache[~mask]==0)), 'stage_n':mask.sum(1).tolist()}
  audit[g]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'file':str(p), 'channels':aa}
  if g=='odd':
   assert np.allclose(cohorts['pairrule']['eve'],dat['eve'],equal_nan=True)
   cohorts['pairrule']['odd']=dat['mrna'];continue
  co='pairrule' if g=='run' else g
  cohorts[co]={'eve':dat['eve'],g:dat['mrna']}
 for co,channels in cohorts.items():
  for ch,a in channels.items():
   for j,s in enumerate(LABELS):
    for k,v in enumerate(a[j]):
     if np.isfinite(v):tidy.append({'cohort':co,'gene':ch,'stage':s,'substage_index':j,'embryo_slot':k+1,'A_minus_B':v})
 audit['shared_pairrule_eve']='Raw eve channels are identical in run/odd workbooks (including NaN); included once. All three channels resampled jointly.'
 pd.DataFrame(tidy).to_csv(ROOT/'temporal_raw_tidy.csv',index=False)
 return cohorts,audit

def boot_cohorts(cohorts,B,suffix=""):
 rng=np.random.default_rng(20261019);curves={};means=[]
 for co,channels in cohorts.items():
  curves[co]={ch:np.zeros((B+1,24)) for ch in channels}
  for j in range(24):
   valid=np.flatnonzero(np.stack([np.isfinite(a[j]) for a in channels.values()]).any(0))
   ix=valid[rng.integers(len(valid),size=(B,len(valid)))]
   for ch,a in channels.items():
    curves[co][ch][0,j]=np.nanmean(a[j])
    with warnings.catch_warnings():
     warnings.simplefilter('ignore',RuntimeWarning);curves[co][ch][1:,j]=np.nanmean(a[j][ix],1)
    vv=a[j,np.isfinite(a[j])];n=len(vv);m=np.mean(vv);sd=np.std(vv,ddof=1)
    means.append(dict(cohort=co,gene=ch,stage=LABELS[j],substage_index=j,n=n,mean=m,sd=sd,sem=sd/np.sqrt(n)))
 pd.DataFrame(means).to_csv(ROOT/f'temporal{suffix}_stage_summary.csv',index=False)
 return curves

def eve_anchors(y):
 # c.1 to c.3, no smoothing. Earliest exact maximum selected; ties recorded.
 anchors=np.empty((len(y),8));ties=np.zeros((len(y),8),bool)
 for c in range(8):
  v=y[:,3*c:3*c+3];valid=np.isfinite(v).all(1);maxv=np.nanmax(v,1)
  anchors[:,c]=np.nanargmax(v,1)+3*c;ties[:,c]=(v==maxv[:,None]).sum(1)>1;anchors[~valid,c]=np.nan
 return anchors,ties

def phase(times,anchors):
 # anchor coordinate is cycle number c. No extrapolation beyond sampled EVE anchors.
 out=np.full(len(times),np.nan)
 for c in range(7):
  a=anchors[:,c];b=anchors[:,c+1];mask=np.isfinite(times)&np.isfinite(a)&np.isfinite(b)&(times>=a)&(times<=b)&(b>a)
  out[mask]=c+1+(times[mask]-a[mask])/(b[mask]-a[mask])
 return out

def gap_metrics(y,g,q,norm):
 lo,hi=[LABELS.index(s) for s in WINDOWS[g]];baseline=np.nanmin(y,1);pk=np.nanmax(y[:,lo:hi+1],1) if norm=='episode' else np.nanmax(y,1);thr=baseline+q*(pk-baseline)
 above=y>=thr[:,None];out={};r=np.arange(len(y));finite=np.isfinite(y[:,lo:hi+1]).all(1)
 for dr in ['rise','fall']:
  cc=((~above[:,lo:hi])&above[:,lo+1:hi+1]) if dr=='rise' else (above[:,lo:hi]&(~above[:,lo+1:hi+1]))
  i=(cc.argmax(1) if dr=='rise' else cc.shape[1]-1-cc[:,::-1].argmax(1))+lo
  with np.errstate(divide='ignore',invalid='ignore'):v=i+(thr-y[r,i])/(y[r,i+1]-y[r,i])
  v[(~cc.any(1))|(~finite)]=np.nan;out[dr]=v
 return out

def pulse_metrics(y,left,right,q):
 # Same common c.1-c.3 peak window for every gene. No allowed ordering is hard-coded.
 r=np.arange(len(y));peak=np.argmax(y[:,left+1:right],1)+left+1;p=y[r,peak]
 finite=np.isfinite(y[:,left:right+1]).all(1)
 # A boundary candidate with a higher adjacent outside-cycle sample is not a resolved pulse maximum.
 resolved=(p>=y[r,peak-1])&(p>=y[r,peak+1])&finite
 localmin=np.zeros_like(y,dtype=bool)
 localmin[:,1:-1]=(y[:,1:-1]<=y[:,:-2])&(y[:,1:-1]<=y[:,2:])&((y[:,1:-1]<y[:,:-2])|(y[:,1:-1]<y[:,2:]))
 xx=np.arange(24)[None,:]
 lc=localmin&(xx<peak[:,None])&(xx>=peak[:,None]-3)
 rc=localmin&(xx>peak[:,None])&(xx<=peak[:,None]+3)
 li=23-lc[:,::-1].argmax(1);ri=rc.argmax(1)
 bracketed=lc.any(1)&rc.any(1)&resolved
 out={'peak':peak.astype(float)};out['peak'][~resolved]=np.nan
 for dr in ['rise','fall']:
  b=y[r,li] if dr=='rise' else y[r,ri];threshold=b+q*(p-b);above=y>=threshold[:,None]
  cc=((~above[:,:-1])&above[:,1:]) if dr=='rise' else (above[:,:-1]&(~above[:,1:]))
  ixs=np.arange(23)[None,:]
  region=(ixs>=li[:,None])&(ixs<peak[:,None]) if dr=='rise' else (ixs>=peak[:,None])&(ixs<ri[:,None])
  cc&=region;i=22-cc[:,::-1].argmax(1) if dr=='rise' else cc.argmax(1)
  with np.errstate(divide='ignore',invalid='ignore'):v=i+(threshold-y[r,i])/(y[r,i+1]-y[r,i])
  v[(~cc.any(1))|(~bracketed)|(p<=b)]=np.nan;out[dr]=v
 return out

def centroid_anchors(y):
 # Exact integrals of the piecewise-linear EVE signal above its cycle-window minimum.
 # First cycle lacks preceding trough, so it is left undefined instead of extrapolated.
 out=np.full((len(y),8),np.nan)
 for c in range(1,8):
  lo=3*c-1;hi=3*c+2;v=y[:,lo:hi+1];w=v-np.nanmin(v,1)[:,None]
  a=w[:,:-1];b=w[:,1:];x=np.arange(lo,hi)[None,:]
  den=np.sum((a+b)/2,1);num=np.sum(x*(a+b)/2+(a+2*b)/6,1)
  with np.errstate(divide='ignore',invalid='ignore'):out[:,c]=num/den
  out[(~np.isfinite(v).all(1))|(den<=0),c]=np.nan
 return out

def addrow(rows,meta,unit,v):
 row={**meta,'unit':unit,**summary(v)}
 if unit=='ordinal_substage' and 'contrast' not in meta:row['bracket_start'],row['bracket_end']=brackets(v[0])
 rows.append(row)

def main(B):
 cohorts,audit=load_sources();curves=boot_cohorts(cohorts,B)
 paired={}
 for co,cc in cohorts.items():
  keep=np.stack([np.isfinite(a) for a in cc.values()]).all(0)
  paired[co]={g:np.where(keep,a,np.nan) for g,a in cc.items()}
  audit.setdefault('complete_pairs',{})[co]=int(keep.sum())
 regcurves=boot_cohorts(paired,B,'_paired');anchors={};centroids={};anchortable=[];anchoraudit={};features={};featuremeta=[]
 def savefeature(key,meta,rawtimes,regtimes,co):
  features[key+'_ordinal']=rawtimes;features[key+'_paired_ordinal']=regtimes
  features[key+'_raw_peak_phase']=phase(regtimes,anchors[co]);features[key+'_centroid_phase']=phase(regtimes,centroids[co])
  featuremeta.append({'key':key,'cohort':co,**meta})

 for co,cc in regcurves.items():
  anc,ties=eve_anchors(cc['eve']);anchors[co]=anc;centroids[co]=centroid_anchors(cc['eve']);anchoraudit[co]={}
  for c in range(8):
   freq={LABELS[i]:float(np.mean(anc[1:,c]==i)) for i in range(c*3,c*3+3)}
   anchoraudit[co][str(c+1)]={'observed_peak_stage':LABELS[int(anc[0,c])],'bootstrap_peak_stage_fractions':freq,'bootstrap_exact_tie_fraction':float(ties[1:,c].mean())}
   anchortable.append(dict(cohort=co,cycle=c+1,observed_peak_stage=LABELS[int(anc[0,c])],exact_tie_fraction=float(ties[1:,c].mean()),**{f'fraction_substage_{k+1}':freq[LABELS[3*c+k]] for k in range(3)}))
 pd.DataFrame(anchortable).to_csv(ROOT/'temporal_eve_anchors.csv',index=False)
 # Global pooled eve, deduped, is display/reference audit only; statistical registration uses same-cohort eve.
 pooled=np.concatenate([v['eve'] for v in cohorts.values()],axis=1)
 poolmean=np.nanmean(pooled,1);poolsem=np.nanstd(pooled,axis=1,ddof=1)/np.sqrt(np.isfinite(pooled).sum(1))
 pd.DataFrame({'stage':LABELS,'substage_index':range(24),'n':np.isfinite(pooled).sum(1),'mean':poolmean,'sem':poolsem}).to_csv(ROOT/'temporal_pooled_eve.csv',index=False)
 gaprows=[];gapcontrast=[];gapstore={};peakrows=[]
 for norm in ['episode','global']:
  for q in [.25,.5,.75]:
   mm={g:gap_metrics(curves[g][g],g,q,norm) for g in GAPS};pm={g:gap_metrics(regcurves[g][g],g,q,norm) for g in GAPS}
   for g,metrics in mm.items():
    for event,times in metrics.items():
     meta=dict(gene=g,event=event,normalization=norm,threshold=q,window_start=WINDOWS[g][0],window_end=WINDOWS[g][1]);addrow(gaprows,{**meta,'anchor_method':'none'},'ordinal_substage',times)
     addrow(gaprows,{**meta,'anchor_method':'raw_peak'},'eve_cycle_phase',phase(pm[g][event],anchors[g]));addrow(gaprows,{**meta,'anchor_method':'signal_centroid'},'eve_cycle_phase',phase(pm[g][event],centroids[g]))
     savefeature(f'gap_{g}_{norm}_{int(q*100)}_{event}',meta,times,pm[g][event],g)
   for a,b in [('KR','SVB'),('MLPT','KR'),('MLPT','SVB'),('GT','MLPT'),('GT','KR'),('GT','SVB')]:
    for event in ['rise','fall']:
     meta=dict(contrast=f'{a}_minus_{b}',event=event,normalization=norm,threshold=q)
     for am,unit,values in [('none','ordinal_substage',mm[a][event]-mm[b][event]),('raw_peak','eve_cycle_phase',phase(pm[a][event],anchors[a])-phase(pm[b][event],anchors[b])),('signal_centroid','eve_cycle_phase',phase(pm[a][event],centroids[a])-phase(pm[b][event],centroids[b]))]:
      addrow(gapcontrast,{**meta,'anchor_method':am},unit,values)
   gapstore[f'{norm}_{int(100*q)}']=mm
 # Actual raw sampled maxima for early episodes; gt treated explicitly as two peaks, not merged peak ordering.
 for g,window,label in [('KR',('1.1','4.3'),'first'),('SVB',('1.1','4.3'),'first'),('MLPT',('1.3','5.1'),'first'),('GT',('3.1','4.1'),'pulse1'),('GT',('4.1','4.3'),'pulse2')]:
  lo,hi=[LABELS.index(s) for s in window];v=curves[g][g][:,lo:hi+1];times=np.argmax(v,1).astype(float)+lo;times[~np.isfinite(v).all(1)]=np.nan
  freq={LABELS[i]:float(np.mean(times[1:]==i)) for i in range(lo,hi+1)}
  pv=regcurves[g][g][:,lo:hi+1];ptimes=np.argmax(pv,1).astype(float)+lo;ptimes[~np.isfinite(pv).all(1)]=np.nan
  meta=dict(gene=g,event=label,observed_peak_stage=LABELS[int(times[0])],bootstrap_stage_fractions=json.dumps(freq))
  for am,unit,vals in [('none','ordinal_substage',times),('raw_peak','eve_cycle_phase',phase(ptimes,anchors[g])),('signal_centroid','eve_cycle_phase',phase(ptimes,centroids[g]))]:addrow(peakrows,{**meta,'anchor_method':am},unit,vals)
  savefeature(f'gap_{g}_{label}_peak',meta,times,ptimes,g)
 pd.DataFrame(gaprows).to_csv(ROOT/'temporal_gap_events.csv',index=False);pd.DataFrame(gapcontrast).to_csv(ROOT/'temporal_gap_contrasts.csv',index=False);pd.DataFrame(peakrows).to_csv(ROOT/'temporal_gap_peaks.csv',index=False)
 # Pairrule cycles2-4; fixed windows refer to sampling substages and are identical across bootstrap draws.
 prrows=[];prcontrast=[];censor=[]
 for cycle in [2,3,4]:
  windows={g:(3*cycle-4,3*cycle) for g in ['eve','run','odd']}
  for q in [.25,.5,.75]:
   mm={g:pulse_metrics(curves['pairrule'][g],*windows[g],q) for g in ['eve','run','odd']}
   for g,metrics in mm.items():
    for event,times in metrics.items():
     if event=='peak' and q!=.5:continue
     obs=curves['pairrule'][g][0];l,r=windows[g];pk=l+1+int(np.argmax(obs[l+1:r])); mins=[i for i in range(1,23) if obs[i]<=obs[i-1] and obs[i]<=obs[i+1] and (obs[i]<obs[i-1] or obs[i]<obs[i+1])]
     leftmins=[i for i in mins if pk-3<=i<pk];rightmins=[i for i in mins if pk<i<=pk+3];lt=max(leftmins) if leftmins else None;rt=min(rightmins) if rightmins else None
     required=max(r-1,pk+1,(rt+1 if event=='fall' and rt is not None else 0))
     meta=dict(gene=g,cycle=cycle,event=event,threshold=q,window_start=LABELS[l+1],window_end=LABELS[r-1],left_trough_stage=LABELS[lt] if lt is not None else '',right_trough_stage=LABELS[rt] if rt is not None else '',required_through_stage=LABELS[required],peak_at_common_window_boundary=bool(pk in [l+1,r-1]),pulse_amplitude_left=float(obs[pk]-obs[lt]) if lt is not None else np.nan,pulse_amplitude_right=float(obs[pk]-obs[rt]) if rt is not None else np.nan)
     addrow(prrows,{**meta,'anchor_method':'none'},'ordinal_substage',times)
     addrow(prrows,{**meta,'anchor_method':'raw_peak'},'eve_cycle_phase',phase(times,anchors['pairrule']));addrow(prrows,{**meta,'anchor_method':'signal_centroid'},'eve_cycle_phase',phase(times,centroids['pairrule']))
     savefeature(f'pairrule_{g}_cycle{cycle}_{int(q*100)}_{event}',meta,times,times,'pairrule')
     for stage in ['3.3','4.1','4.3','5.2']:
      end=LABELS.index(stage);valid=np.isfinite(times[1:]);censor.append({**meta,'comparison_stage':stage,'observed_event_completed':bool(np.isfinite(times[0]) and times[0]<=end),'bootstrap_event_completed_fraction':float(np.mean(np.isfinite(times[1:])&(times[1:]<=end))),'event_extraction_supported_by_this_stage':bool(required<=end),'note':'Event may occur before comparison stage but its peak/trough classification or amplitude estimate can require later observations.'})
   for a,b in [('run','eve'),('odd','run'),('odd','eve')]:
    for event in ['rise','peak','fall']:
     if event=='peak' and q!=.5:continue
     for am,unit,v in [('none','ordinal_substage',mm[a][event]-mm[b][event]),('raw_peak','eve_cycle_phase',phase(mm[a][event],anchors['pairrule'])-phase(mm[b][event],anchors['pairrule'])),('signal_centroid','eve_cycle_phase',phase(mm[a][event],centroids['pairrule'])-phase(mm[b][event],centroids['pairrule']))]:addrow(prcontrast,dict(contrast=f'{a}_minus_{b}',cycle=cycle,event=event,threshold=q,anchor_method=am),unit,v)
 pd.DataFrame(prrows).to_csv(ROOT/'temporal_pairrule_events.csv',index=False);pd.DataFrame(prcontrast).to_csv(ROOT/'temporal_pairrule_contrasts.csv',index=False);pd.DataFrame(censor).to_csv(ROOT/'temporal_pairrule_stage_availability.csv',index=False)
 # Tipping offsets: deterministic sensitivity, no invented probability distribution or jitter CI.
 tip=[]
 for row in gapcontrast:
  if row['threshold']!=.5:continue
  est=row['estimate'];lo=row['ci95_low'];hi=row['ci95_high']
  if est is None:continue
  bound=min(abs(lo),abs(hi)) if lo*hi>0 else 0.
  tip.append({k:row[k] for k in ['contrast','event','normalization','threshold','unit','anchor_method']}|{'signed_residual_offset_to_null_observed':-est,'absolute_residual_offset_to_null_observed':abs(est),'absolute_adverse_offset_to_make_ci_include_zero':bound,'note':'Subtract this contrast-specific systematic residual offset; not an estimated staging error. EVE phase cancels a uniform common shift of gene and reference within each cohort, but does not remove within-stage or nonlinear timing uncertainty.'})
 pd.DataFrame(tip).to_csv(ROOT/'temporal_tipping_offsets.csv',index=False)
 # Source-backed figures: independent fixed embryos, stage mean and SEM, no fitted/smoothed curves.
 for name,genes,cohortmap in [('gap',GAPS,{g:g for g in GAPS}),('pairrule',['eve','run','odd'],{g:'pairrule' for g in ['eve','run','odd']})]:
  fig,axes=plt.subplots(len(genes),1,figsize=(11,2.25*len(genes)),sharex=True,layout='constrained')
  for ax,g in zip(axes,genes):
   a=cohorts[cohortmap[g]][g];mean=np.nanmean(a,1);sem=np.nanstd(a,1,ddof=1)/np.sqrt(np.isfinite(a).sum(1));x=np.arange(24)
   ax.errorbar(x,mean,yerr=sem,fmt='o-',c=COLORS[g],lw=1.5,ms=4,capsize=2,label='Stage mean ± SEM')
   for i in x:
    vals=a[i,np.isfinite(a[i])];ax.scatter(i+np.linspace(-.13,.13,len(vals)),vals,s=12,color=COLORS[g],alpha=.38)
   ax.set_ylabel(f'{g}\nA − B intensity');ax.spines[['top','right']].set_visible(False);ax.set_xlim(-.4,23.4)
   for k in range(1,8):ax.axvline(3*k-.5,c='#dddddd',lw=.7)
  axes[0].set_title('Unsmoothed posterior temporal measurements: '+('gap genes' if name=='gap' else 'eve, run and odd'),loc='left',weight='bold');axes[0].legend(loc='upper left',fontsize=9)
  axes[-1].set_xticks(range(24),LABELS,rotation=45);axes[-1].set_xlabel('Ordered developmental substages (not calibrated minutes)');fig.savefig(ROOT/f'temporal_{name}_raw.png',dpi=180);plt.close(fig)
 fig,ax=plt.subplots(figsize=(11,4),layout='constrained')
 for co,cc in curves.items():
  y=cc['eve'][0];z=(y-y.min())/(y.max()-y.min());ax.plot(range(24),z,'o-',ms=3,lw=1,alpha=.75,label=co)
 ax.set_xticks(range(24),LABELS,rotation=45);ax.set_ylabel('EVE stage mean / cohort-wide range');ax.set_xlabel('Ordered developmental substage');ax.set_title('EVE registration audit: independently sampled cohorts',loc='left',weight='bold');ax.legend(ncol=6,fontsize=8);ax.spines[['top','right']].set_visible(False);fig.savefig(ROOT/'temporal_eve_cohort_audit.png',dpi=180);plt.close(fig)
 methods={'bootstrap_replicates':B,'seed':20261019,'source_audit':audit,'eve_anchor_audit':anchoraudit,'scope':'Four gap genes and pairrule eve/run/odd. HB contributes only its independently measured eve channel to reference audit/pool.', 'replication':'Columns are embryo slots at each stage, not trajectories. All-observed ordinal analyses include all observed target measurements. Registered analyses use complete multichannel pairs; KR115 pairs versus116 target values. Pair-rule106 complete triples. Resample whole within-stage records, distinct gap cohorts independently. No zero imputation.', 'smoothing':'None. Means are joined linearly for display/interpolated thresholds. Sampled maxima remain coarse sampled-stage events.', 'gap_events':{'windows':WINDOWS,'baseline':'whole-series minimum stage mean','thresholds':[.25,.5,.75],'normalizations':['first episode maximum','global maximum'],'gt':'two pulses treated as envelope for rise/fall; maxima assessed separately','rise':'first upward crossing in window','fall':'last downward crossing in window'},'pairrule_events':'Cycles2-4: common c.1-c.3 maximum window for ALL three genes allows any order. Candidate must not be lower than adjacent outside-window sample. Rise/fall use nearest raw local minima on each side within one cycle; undefined if absent. Each side has its own trough-to-peak threshold. Narrow or zero CIs can be structurally induced by interpolation between adjacent samples and do not quantify within-substage timing.', 'eve_registration':'Cohort-specific raw eve maximum within c.1-c.3 is anchor c. Piecewise-linear fractional cycle coordinate between successive anchors, no extrapolation. Exact ties choose earliest sampled maximum and tie frequency is reported. This is one explicitly uncertain registration, not a latent-time reconstruction.', 'uncertainty':'Percentile95% bootstrap intervals conditional on assigned stages, windows, threshold definitions and channel pairing. Not simultaneous across all comparisons. Does not estimate stage classification error or nonlinear stage durations.', 'offset_sensitivity':'Report residual differential offsets needed to null observed contrasts or make CI include zero. No assumed timing jitter was injected. Uniform cohort shift shared by gene+eve cancels in relative phase if both landmarks shift together.', 'completion':'Retrospective pulse definitions can require stages later than comparison stage; availability table flags these rather than pretending future data existed at that stage.'}
 # Feature vectors allow external spatial comparisons using the exact same bootstrap draws.
 np.savez_compressed(ROOT/'temporal_bootstrap_features.npz',**features)
 (ROOT/'temporal_bootstrap_feature_index.json').write_text(json.dumps(featuremeta,indent=2))
 peak_audit=[]
 for fm in featuremeta:
  if not fm['key'].startswith('pairrule_') or fm['event']!='peak':continue
  tt=features[fm['key']+'_ordinal'];vals,counts=np.unique(tt[1:][np.isfinite(tt[1:])],return_counts=True)
  peak_audit.append({**fm,'observed_peak_stage':LABELS[int(tt[0])] if np.isfinite(tt[0]) else '', 'valid_bootstrap_fraction':float(np.isfinite(tt[1:]).mean()),'bootstrap_sampled_peak_frequencies':json.dumps({LABELS[int(v)]:int(n)/B for v,n in zip(vals,counts)}),'interpretation':'A shared sampled maximum does not establish simultaneous true maxima. Boundary/window and unresolved fraction must be retained.'})
 pd.DataFrame(peak_audit).to_csv(ROOT/'temporal_pairrule_peak_audit.csv',index=False)
 # Ordinal event relationships to sampled EVE landmarks. Centroids are a coordinate-dependent sensitivity, not invariant to nonlinear time warps.
 relation=[]
 for meta in featuremeta:
  key=meta['key'];co=meta['cohort'];t=features[key+'_paired_ordinal']
  for am,aa in [('raw_peak',anchors[co]),('signal_centroid',centroids[co])]:
   for c in range(8):
    d=t-aa[:,c];valid=np.isfinite(d[1:]);vv=d[1:][valid]
    if not np.isfinite(d[0]):continue
    relation.append({**meta,'anchor_method':am,'eve_cycle':c+1,**summary(d),'bootstrap_fraction_before':float(np.mean(vv<0)) if len(vv) else np.nan,'bootstrap_fraction_equal_sampled_coordinate':float(np.mean(vv==0)) if len(vv) else np.nan,'bootstrap_fraction_after':float(np.mean(vv>0)) if len(vv) else np.nan})
 pd.DataFrame(relation).to_csv(ROOT/'temporal_event_eve_relations.csv',index=False)
 cr=[]
 for co,aa in centroids.items():
  for c in range(8):cr.append({'cohort':co,'cycle':c+1,**summary(aa[:,c])})
 pd.DataFrame(cr).to_csv(ROOT/'temporal_eve_centroids.csv',index=False)
 methods['eve_centroids']='Sensitivity anchors: exact integral of piecewise-linear raw EVE signal above local minimum in previous cycle .3 through current cycle .3. First cycle undefined because no preceding trough. Fractional distances remain conditional on assumed within-cycle coordinate map.'
 (ROOT/'temporal_methods_and_audit.json').write_text(json.dumps(methods,indent=2))
 np.savez_compressed(ROOT/'temporal_bootstrap_curves.npz',**{f'{co}_{ch}':v for co,cc in curves.items() for ch,v in cc.items()},**{f'{co}_eve_anchors':a for co,a in anchors.items()},**{f'{co}_eve_centroids':a for co,a in centroids.items()},**{f'paired_{co}_{ch}':v for co,cc in regcurves.items() for ch,v in cc.items()})
 print('Completed temporal direct analysis; B=',B)
 print(pd.DataFrame(gapcontrast).query("normalization=='episode' and threshold==0.5")[['contrast','event','unit','estimate','ci95_low','ci95_high']].to_string(index=False))
 print(pd.DataFrame(prcontrast).query('threshold==0.5')[['contrast','cycle','event','unit','estimate','ci95_low','ci95_high']].to_string(index=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--bootstraps',type=int,default=20000);main(p.parse_args().bootstraps)
