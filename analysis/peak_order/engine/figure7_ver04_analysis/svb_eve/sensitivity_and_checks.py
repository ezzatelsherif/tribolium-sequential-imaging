from pathlib import Path
import sys,json
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));from extract_svb_eve import smoothing,v03,DATA
raw=pd.read_csv(DATA/'spatial_raw_profiles.csv');eve=pd.read_csv(HERE/'spatial_features.csv');eve=eve[eve.gene=='eve'].set_index(['stage','embryo'])
rows=[]
for (stage,emb),p in raw[(raw.gene=='svb')&raw.embryo.isin([1,2,3])].groupby(['stage','embryo']):
 p=p.sort_values('x_ap');x=p.x_ap.to_numpy();y=p.raw_intensity.to_numpy()
 for frac in [.05,.1,.15,.2]:
  ys=smoothing.smooth_signal(y,frac=frac,polyorder=3);f,peaks,cross=v03.spatial_feature(x,ys,.2,1)
  for event in ['rise','peak','fall']:
   e=float(eve.loc[(stage,emb),event]);v=float(f[event]);rows.append(dict(stage=stage,embryo=emb,smoothing_fraction=frac,event=event,svb=v,eve=e,svb_minus_eve=v-e,svb_peak_at_window_edge=f['peak_at_window_edge'],svb_multipeak=f['multi_peak_flag'],excluded_neural=stage=='Stage5_2'))
s=pd.DataFrame(rows);s.to_csv(HERE/'spatial_svb_smoothing_sensitivity.csv',index=False)
r=[]
for (stage,event),p in s.groupby(['stage','event']):
 for frac,pp in p.groupby('smoothing_fraction'):
  r.append(dict(stage=stage,event=event,fraction=frac,negative=int((pp.svb_minus_eve<0).sum()),positive=int((pp.svb_minus_eve>0).sum()),missing=int(pp.svb_minus_eve.isna().sum())))
r=pd.DataFrame(r);r.to_csv(HERE/'spatial_svb_smoothing_sensitivity_summary.csv',index=False);print(r.to_string(index=False))
# Verify affine normalization does not change features under the exact original smoother.
checks=[]
for (stage,emb,gene),p in raw[(raw.gene.isin(['svb','kr','eve']))&raw.embryo.isin([1,2,3])].groupby(['stage','embryo','gene']):
 p=p.sort_values('x_ap');x=p.x_ap.to_numpy();y=p.raw_intensity.to_numpy()
 if gene=='svb':
  normalized=smoothing.smooth_signal(smoothing.normalize_range(y),frac=.1,polyorder=3);unscaled=smoothing.smooth_signal(y,frac=.1,polyorder=3)
 else:
  _,normalized=smoothing.spatial(x,y,stage,emb,gene);_,unscaled=smoothing.spatial(x,y,stage,emb,gene,normalize=False)
 err=np.max(np.abs(normalized-(unscaled-y.min())/np.ptp(y)));assert err<1e-10
 checks.append(dict(stage=stage,embryo=int(emb),gene=gene,maximum_affine_smoothing_error=float(err)))
json.dump({'normalization_commutation_checks':checks,'assertion':'SG smoothing raw intensity is affine-equivalent to original normalization then SG; half-height normalization is applied after smoothing.'},open(HERE/'validation.json','w'),indent=2)
