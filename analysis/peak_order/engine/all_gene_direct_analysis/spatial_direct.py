"""Unsmoothed spatial audit of gap and pair-rule genes, original n=3/stage.

Curated coordinates identify coarse expression episodes only. All reported
measurements are recomputed from raw CSV samples, with piecewise-linear
interpolation for thresholds/integration. No Savitzky-Golay or other smoothing.
"""
from pathlib import Path
import json, itertools
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import t,binomtest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
SRC=ROOT
CFGS=json.loads((SRC/'reconstructed_points.json').read_text())['configs']
GENES=['kr','svb','mlpt','gt','eve','runt','odd']
COLORS={'kr':'#c93631','svb':'#228bd0','mlpt':'#258b48','gt':'#bd8b15','eve':'#333333','runt':'#a13b96','odd':'#ec7c25'}

def cut(x,y,lo,hi):
    z=(x>lo)&(x<hi)
    return np.r_[lo,x[z],hi],np.r_[np.interp(lo,x,y),y[z],np.interp(hi,x,y)]

def crossing(x,y,q,side):
    idx=np.flatnonzero(y>=q)
    if not len(idx): return np.nan,True
    k=int(idx[0] if side=='rise' else idx[-1])
    if (side=='rise' and k==0) or (side=='fall' and k==len(x)-1): return np.nan,True
    i,j=(k-1,k) if side=='rise' else (k,k+1)
    return float(x[i]+(q-y[i])*(x[j]-x[i])/(y[j]-y[i])),False

def centroid(x,y):
    # Exact integral of x*y for each piecewise-linear segment after local min.
    y=y-y.min(); dx=np.diff(x)
    area=np.sum(dx*(y[:-1]+y[1:])/2)
    moment=np.sum(dx*((2*x[:-1]+x[1:])*y[:-1]+(x[:-1]+2*x[1:])*y[1:])/6)
    return float(moment/area) if area>0 else np.nan

def episode(x,y,lo,hi):
    xx,yy=cut(x,y,lo,hi)
    # Select maxima only from original observed samples, never a resampling grid.
    ids=np.flatnonzero((x>=lo)&(x<=hi)); k=int(ids[np.argmax(y[ids])])
    xp=float(x[k]); yp=float(y[k]); loc=int(np.searchsorted(xx,xp))
    left_min=float(yy[:loc+1].min());right_min=float(yy[loc:].min())
    peaks,props=find_peaks(yy,prominence=.15*max(float(np.ptp(yy)),1e-12))
    prominent=[float(xx[i]) for i in peaks]
    split=any(abs(a-b)>=.04 for a,b in itertools.combinations(prominent,2))
    ans={'peak':xp,'centroid':centroid(xx,yy),'window_lo':lo,'window_hi':hi,
         'boundary_definition':'coarse_episode_envelope',
         'peak_height':yp,'left_min':left_min,'right_min':right_min,
         'left_min_at_window_edge':bool(np.argmin(yy[:loc+1])==0),
         'right_min_at_window_edge':bool(np.argmin(yy[loc:])+loc==len(yy)-1),
         'multi_peak_flag':split,'prominent_peak_positions':json.dumps(prominent)}
    ans['peak_at_window_edge']=bool(k==ids[0] or k==ids[-1])
    for q in [.25,.5,.75]:
        # Side-specific endpoint baseline, matching temporal direct analysis.
        level_l=left_min+q*(yp-left_min);level_r=right_min+q*(yp-right_min)
        # Boundary of entire episode on each side, retaining all above-level lobes.
        a,ca=crossing(xx[:loc+1],yy[:loc+1],level_l,'rise')
        b,cb=crossing(xx[loc:],yy[loc:],level_r,'fall')
        ans[f'rise_{q}']=a;ans[f'fall_{q}']=b
        ans[f'rise_{q}_censored']=ca;ans[f'fall_{q}_censored']=cb
        ans[f'rise_{q}_level']=level_l;ans[f'fall_{q}_level']=level_r
        # Connected support around the selected peak, using the SAME coarse
        # episode baselines. Pixel-scale local minima are not new pulse baselines.
        left_below=np.flatnonzero(yy[:loc]<level_l)
        right_below=np.flatnonzero(yy[loc+1:]<level_r)+loc+1
        if len(left_below):
            k0=int(left_below[-1]);ans[f'connected_rise_{q}']=float(xx[k0]+(level_l-yy[k0])*(xx[k0+1]-xx[k0])/(yy[k0+1]-yy[k0]))
        else:ans[f'connected_rise_{q}']=np.nan
        if len(right_below):
            k1=int(right_below[0]);ans[f'connected_fall_{q}']=float(xx[k1-1]+(level_r-yy[k1-1])*(xx[k1]-xx[k1-1])/(yy[k1]-yy[k1-1]))
        else:ans[f'connected_fall_{q}']=np.nan
    return ans

def phase(x,anchors):
    aa=sorted((int(c),float(v)) for c,v in anchors.items() if np.isfinite(v))
    if not aa or not np.isfinite(x):return np.nan,'unavailable'
    if x<aa[0][1]:return np.nan,f'before_{aa[0][0]}'
    if x>aa[-1][1]:return np.nan,f'after_{aa[-1][0]}'
    for (c,a),(d,b) in zip(aa[:-1],aa[1:]):
        if a<=x<=b:return c+(x-a)/(b-a)*(d-c),f'{c}_to_{d}'
    return float(aa[0][0]),f'at_{aa[0][0]}'

def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    profiles={};raw_rows=[];domain=[];events=[];anchors=[]
    for cfg in CFGS:
        stage=cfg['stage']
        for emb in [1,2,3]:
            d=pd.read_csv(SRC/'data'/stage/f'EMB_{emb}.csv')
            xx=d.iloc[:,0].to_numpy(float);x=(xx-xx.min())/np.ptp(xx)
            profiles[(stage,emb)]={}
            for g in GENES:
                raw=d[g].to_numpy(float);y=(raw-raw.min())/np.ptp(raw)
                profiles[(stage,emb)][g]=(x,y)
                raw_rows.append(pd.DataFrame({'stage':stage,'embryo':emb,'gene':g,'x_ap':x,
                    'distance_from_posterior_normalized':1-x,'raw_intensity':raw,'normalized_intensity':y}))
                lo=.45 if g in ['gt','mlpt'] else .2
                ids=np.flatnonzero(x>=lo);kp=int(ids[np.argmax(y[ids])])
                for q in [.25,.5,.75]:
                    xc,yc=cut(x,y,lo,1.)
                    rise,cr=crossing(xc,yc,q,'rise');fall,cf=crossing(xc,yc,q,'fall')
                    domain.append({'stage':stage,'embryo':emb,'gene':g,'threshold_fraction':q,
                        'window_lo':lo,'window_hi':1.,'peak':float(x[kp]),'rise':rise,'fall':fall,
                        'rise_censored':cr,'fall_censored':cf,'excluded_neural':g in ['svb','kr'] and stage=='Stage5_2'})
            # Pair-rule cycles: even historical landmarks identify peaks; odd
            # landmarks bound episodes. Their numerical positions are not outcomes.
            for g in ['eve','runt','odd']:
                x,y=profiles[(stage,emb)][g]
                pts=cfg['points'][str(emb)][g]
                for i in range(1,len(pts),2):
                    cyc=i//2+2
                    lo=float(pts[i-1]);has_right=i+1<len(pts)
                    hi=float(pts[i+1]) if has_right else 1.
                    met=episode(x,y,lo,hi)
                    rec={'stage':stage,'embryo':emb,'gene':g,'group':'pair_rule',
                        'event':f'{g}_cycle{cyc}','cycle':cyc,'complete_coarse_window':has_right,
                        'excluded_neural':False,**met}
                    events.append(rec)
                    if g=='eve':
                        anchors.append({**rec,'posterior_developing_cycle':stage=='Stage4_1' and cyc==4,
                            'peak_anchor_primary_eligible':not met['multi_peak_flag'] and not(stage=='Stage4_1' and cyc==4)})
            # Gap main episodes and separate gt lobes. svb has no historical
            # landmarks, so its broad posterior domain uses x=.2..1.
            for g in ['kr','mlpt','gt','svb']:
                x,y=profiles[(stage,emb)][g]
                if g=='svb':windows=[('svb_main',.2,1.,False)]
                else:
                    pts=cfg['points'][str(emb)][g]
                    windows=[]
                    for i in range(1,len(pts),2):
                        lo=float(pts[i-1]);complete=i+1<len(pts)
                        hi=float(pts[i+1]) if complete else 1.
                        windows.append((f'{g}_'+('main' if g!='gt' else f'peak{i//2+1}'),lo,hi,complete))
                for name,lo,hi,complete in windows:
                    events.append({'stage':stage,'embryo':emb,'gene':g,'group':'gap',
                        'event':name,'cycle':np.nan,'complete_coarse_window':complete,
                        'excluded_neural':g in ['svb','kr'] and stage=='Stage5_2',**episode(x,y,lo,hi)})
    pd.concat(raw_rows,ignore_index=True).to_csv(ROOT/'spatial_raw_profiles.csv',index=False)
    dd=pd.DataFrame(domain);dd.to_csv(ROOT/'spatial_domain_metrics.csv',index=False)
    ee=pd.DataFrame(events);ee.to_csv(ROOT/'spatial_episode_metrics.csv',index=False)
    aa=pd.DataFrame(anchors);aa.to_csv(ROOT/'spatial_eve_anchors.csv',index=False)
    lobe_control=[]
    for emb in [1,2,3]:
        ar=aa[(aa.stage=='Stage4_3')&(aa.embryo==emb)&(aa.cycle==2)].iloc[0]
        lobes=json.loads(ar.prominent_peak_positions)
        glob=dd[(dd.stage=='Stage4_3')&(dd.embryo==emb)&(dd.gene=='kr')&(dd.threshold_fraction==.5)].iloc[0]['rise']
        loc=ee[(ee.stage=='Stage4_3')&(ee.embryo==emb)&(ee.event=='kr_main')].iloc[0]['rise_0.5']
        for norm,v in [('global_profile_range',glob),('side_specific_local',loc)]:
            lobe_control.append({'stage':'Stage4_3','embryo':emb,'normalization':norm,
                'kr_halfheight_rise':v,'eve2_first_prominent_lobe':min(lobes),
                'eve2_last_prominent_lobe':max(lobes),'eve2_centroid':ar.centroid,
                'eve2_dominant_peak':ar.peak,'offset_to_earliest_lobe_ap':v-min(lobes),
                'before_every_prominent_lobe':bool(v<min(lobes)),
                'between_prominent_lobes':bool(min(lobes)<=v<=max(lobes)),
                'interpretation':'Before centroid/dominant maximum does not establish an anchor-independent reversal.'})
    pd.DataFrame(lobe_control).to_csv(ROOT/'spatial_kr_eve2_all_lobes_control.csv',index=False)
    # Gene-blind window sensitivity: every pair-rule gene searched in the same
    # EVE coarse-cycle interval. Clipped peaks cannot support ordering claims.
    common=[]
    for a in anchors:
        for g in ['eve','runt','odd']:
            x,y=profiles[(a['stage'],a['embryo'])][g]
            identified=any(r['stage']==a['stage'] and r['embryo']==a['embryo']
                and r['gene']==g and r['cycle']==a['cycle'] for r in events)
            common.append({'stage':a['stage'],'embryo':a['embryo'],'gene':g,'group':'pair_rule',
                'event':f'{g}_cycle{int(a["cycle"])}','cycle':int(a['cycle']),
                'search_window':'common_eve_cycle','complete_coarse_window':a['complete_coarse_window'],
                'homologous_event_previously_identified':identified,
                'posterior_developing_cycle':a['posterior_developing_cycle'],
                **episode(x,y,a['window_lo'],a['window_hi'])})
    pd.DataFrame(common).to_csv(ROOT/'spatial_common_window_metrics.csv',index=False)
    offsets=[]
    for (stage,emb,q),sub in dd.groupby(['stage','embryo','threshold_fraction']):
        sub=sub.set_index('gene')
        for g,h in [('gt','mlpt'),('svb','kr')]:
            for feature in ['peak','rise','fall']:
                a,b=float(sub.loc[g,feature]),float(sub.loc[h,feature])
                offsets.append({'stage':stage,'embryo':emb,'gene':g,'reference':h,
                    'feature':feature,'threshold_fraction':q,'offset_ap':a-b,
                    'excluded_neural':(g in ['svb','kr'] or h in ['svb','kr']) and stage=='Stage5_2',
                    'either_censored':not(np.isfinite(a) and np.isfinite(b))})
    oo=pd.DataFrame(offsets);oo.to_csv(ROOT/'spatial_direct_offsets.csv',index=False)
    ss=[]
    for key,sub in oo[~oo.excluded_neural].groupby(['stage','gene','reference','feature','threshold_fraction']):
        vals=sub.offset_ap.dropna().to_numpy();n=len(vals)
        sd=float(vals.std(ddof=1)) if n>1 else np.nan
        se=sd/np.sqrt(n) if n else np.nan
        mean=float(vals.mean()) if n else np.nan
        margin=float(t.ppf(.975,n-1)*se) if n>1 else np.nan
        nz=vals[vals!=0]
        ss.append(dict(zip(['stage','gene','reference','feature','threshold_fraction'],key),
                       n_embryos=n,mean_offset_ap=mean,sd=sd,sem=se,
                       mean_t95_low=mean-margin,mean_t95_high=mean+margin,
                       exact_two_sided_sign_p=float(binomtest(sum(nz>0),len(nz)).pvalue) if len(nz) else np.nan))
    pd.DataFrame(ss).to_csv(ROOT/'spatial_direct_offset_summary.csv',index=False)
    # Extra files are not silently pooled with the notebook's original n=3.
    extra=[];extra_offsets=[]
    for cfg in CFGS:
        stage=cfg['stage'];path=SRC/'data'/stage/'EMB_4.csv'
        if not path.exists():continue
        data=pd.read_csv(path);x=data.iloc[:,0].to_numpy(float);x=(x-x.min())/np.ptp(x)
        for g in ['kr','svb','mlpt','gt']:
            y=data[g].to_numpy(float);y=(y-y.min())/np.ptp(y);lo=.45 if g in ['gt','mlpt'] else .2
            ids=np.flatnonzero(x>=lo);xp=float(x[ids[np.argmax(y[ids])]])
            for q in [.25,.5,.75]:
                xc,yc=cut(x,y,lo,1.);rise,cr=crossing(xc,yc,q,'rise');fall,cf=crossing(xc,yc,q,'fall')
                extra.append({'stage':stage,'embryo':4,'gene':g,'threshold_fraction':q,'window_lo':lo,
                    'peak':xp,'rise':rise,'fall':fall,'rise_censored':cr,'fall_censored':cf,
                    'provenance_unverified':True,'included_in_primary_n':False})
    if extra:
        ex=pd.DataFrame(extra);ex.to_csv(ROOT/'spatial_extra_embryo_domain_metrics.csv',index=False)
        for (stage,q),sub in ex.groupby(['stage','threshold_fraction']):
            sub=sub.set_index('gene')
            for g,h in [('gt','mlpt'),('svb','kr')]:
                for feature in ['peak','rise','fall']:
                    v=float(sub.loc[g,feature]-sub.loc[h,feature])
                    extra_offsets.append({'stage':stage,'embryo':4,'gene':g,'reference':h,'feature':feature,
                        'threshold_fraction':q,'offset_ap':v,'provenance_unverified':True,'included_in_primary_n':False})
        pd.DataFrame(extra_offsets).to_csv(ROOT/'spatial_extra_embryo_offsets.csv',index=False)
    phase_rows=[]
    for (stage,emb),sub in ee.groupby(['stage','embryo']):
        av=aa[(aa.stage==stage)&(aa.embryo==emb)&~aa.posterior_developing_cycle]
        for rep in ['peak','centroid']:
            anc=dict(zip(av.cycle,av[rep]))
            def anchor_status(bracket):
                cycles=[int(v) for v in bracket.split('_') if v.isdigit()]
                used=av[av.cycle.isin(cycles)]
                return {'anchor_has_multi_peak':bool(used.multi_peak_flag.any()),
                    'anchor_primary_eligible':bool(len(used) and (rep=='centroid' or used.peak_anchor_primary_eligible.all()))}
            for r in sub.to_dict('records'):
                for feat in ['peak','centroid','rise_0.25','fall_0.25','rise_0.5','fall_0.5','rise_0.75','fall_0.75']:
                    val,bracket=phase(r[feat],anc)
                    phase_rows.append({'stage':stage,'embryo':emb,'gene':r['gene'],'group':r['group'],
                        'event':r['event'],'feature':feat,'x_ap':r[feat],'anchor_representation':rep,
                        'anchor_method':'raw_peak' if rep=='peak' else 'signal_centroid',
                        'normalization':'side_specific_local',
                        **anchor_status(bracket),
                        'eve_phase':val,'eve_bracket':bracket,'multi_peak_flag':r['multi_peak_flag'],
                        'complete_coarse_window':r['complete_coarse_window'],
                        'excluded_neural':r['excluded_neural'],
                        'anchor_cycles':json.dumps(sorted(anc))})
            # Global profile-range boundaries are a separate sensitivity, never
            # silently equated with local episode-relative boundary definitions.
            gd=dd[(dd.stage==stage)&(dd.embryo==emb)&dd.gene.isin(['kr','svb','mlpt','gt'])]
            for r in gd.to_dict('records'):
                for feature in ['rise','fall']:
                    feat=f'{feature}_{r["threshold_fraction"]}'
                    val,bracket=phase(r[feature],anc)
                    phase_rows.append({'stage':stage,'embryo':emb,'gene':r['gene'],'group':'gap',
                        'event':r['gene']+'_whole_domain','feature':feat,'x_ap':r[feature],
                        'anchor_representation':rep,'anchor_method':'raw_peak' if rep=='peak' else 'signal_centroid',
                        'normalization':'global_profile_range','eve_phase':val,'eve_bracket':bracket,
                        **anchor_status(bracket),
                        'multi_peak_flag':np.nan,'complete_coarse_window':not(r['rise_censored'] or r['fall_censored']),
                        'excluded_neural':r['excluded_neural'],'anchor_cycles':json.dumps(sorted(anc))})
    pd.DataFrame(phase_rows).to_csv(ROOT/'spatial_eve_relative_events.csv',index=False)
    # Shared AP grid is only for display of means/SEM. It never changes event values.
    meanrows=[];grid=np.linspace(0,1,1001)
    for cfg in CFGS:
        stage=cfg['stage']
        for g in GENES:
            ys=np.array([np.interp(grid,*profiles[(stage,e)][g]) for e in [1,2,3]])
            meanrows.append(pd.DataFrame({'stage':stage,'gene':g,'x_ap':grid,'mean':ys.mean(0),
                                         'sem':ys.std(0,ddof=1)/np.sqrt(3),'n_embryos':3}))
    mm=pd.concat(meanrows,ignore_index=True);mm.to_csv(ROOT/'spatial_mean_sem.csv',index=False)
    for name,genes in [('gap',['kr','svb','mlpt','gt']),('pair_rule',['eve','runt','odd'])]:
        fig,axs=plt.subplots(4,4,figsize=(16,12),sharex=True,sharey=True,layout='constrained')
        for i,cfg in enumerate(CFGS):
            stage=cfg['stage']
            for j,emb in enumerate([1,2,3]):
                ax=axs[i,j]
                for g in genes:
                    x,y=profiles[(stage,emb)][g]
                    ax.plot(x,y,color=COLORS[g],lw=.9,label=g,
                            alpha=.4 if g in ['svb','kr'] and stage=='Stage5_2' else 1)
                ax.set_title(f'{stage.replace("Stage","").replace("_",".")}, embryo {emb}',fontsize=10)
            ax=axs[i,3]
            for g in genes:
                m=mm[(mm.stage==stage)&(mm.gene==g)]
                ax.fill_between(m.x_ap,m['mean']-m['sem'],m['mean']+m['sem'],color=COLORS[g],alpha=.15)
                ax.plot(m.x_ap,m['mean'],color=COLORS[g],lw=1.2,label=g)
            ax.set_title('Mean ± SEM, n=3',fontsize=10)
        for ax in axs.flat:
            ax.set(xlim=(0,1),ylim=(-.03,1.07));ax.spines[['top','right']].set_visible(False)
        axs[0,0].legend(ncol=len(genes),fontsize=8,loc='upper left',frameon=False)
        fig.supxlabel('Normalized anterior-to-posterior position (posterior tip = 1)')
        fig.supylabel('Raw intensity / each gene’s own whole-profile range')
        fig.suptitle(f'Unsmoothed {name.replace("_","-")} profiles; original three embryos per stage'+
                     ('\nKr and svb at stage 5.2 are displayed for provenance and excluded from primary inference' if name=='gap' else ''),fontsize=13)
        fig.savefig(ROOT/f'spatial_{name}_raw_profiles.png',dpi=170);plt.close(fig)
    # EVE audit: raw curve, coarse window boundaries, raw maxima, centroids.
    fig,axs=plt.subplots(4,3,figsize=(13,10),sharex=True,sharey=True,layout='constrained')
    for i,cfg in enumerate(CFGS):
        for j,emb in enumerate([1,2,3]):
            ax=axs[i,j];stage=cfg['stage'];x,y=profiles[(stage,emb)]['eve'];ax.plot(x,y,color='k',lw=.8)
            sub=aa[(aa.stage==stage)&(aa.embryo==emb)]
            for r in sub.to_dict('records'):
                ax.axvspan(r['window_lo'],r['window_hi'],alpha=.08,color=['#4075b0','#d87921'][int(r['cycle'])%2])
                ax.plot(r['peak'],np.interp(r['peak'],x,y),'o',color='#c93631',ms=4)
                ax.axvline(r['centroid'],color='#228bd0',ls=':',lw=1)
                ax.text(r['centroid'],1.04,f'C{int(r["cycle"])}'+('*' if r['multi_peak_flag'] else ''),ha='center',fontsize=8)
            ax.set_title(f'{stage.replace("Stage","").replace("_",".")}, embryo {emb}',fontsize=10)
            ax.set(ylim=(-.02,1.14),xlim=(0,1));ax.spines[['top','right']].set_visible(False)
    fig.suptitle('EVE anchor audit: red = raw maximum, blue dotted = raw-signal centroid\n* multiple prominent raw peaks; shaded windows identify historical cycles only',fontsize=12)
    fig.supxlabel('Normalized AP position');fig.supylabel('Normalized raw EVE intensity')
    fig.savefig(ROOT/'spatial_eve_anchor_audit.png',dpi=170);plt.close(fig)
    (ROOT/'spatial_methods.json').write_text(json.dumps({
        'primary_embryos':[1,2,3],'stages':[c['stage'] for c in CFGS],
        'genes':GENES,'normalization':'Each gene/embryo whole-profile min-max; AP position=(distance-min)/(max-min).',
        'smoothing':'None. Linear interpolation only for crossings and display grid; exact linear-segment integration for centroids.',
        'raw_peaks':'Actual observed CSV samples, not interpolated grid maxima. Exact ties use earliest sample.',
        'display_means':'Unshifted normalized AP grid; biological n=3. Interpolation is for display only.',
        'global_domain_metrics':'Kr/svb x=.2..1; gt/mlpt x=.45..1 to exclude separate anterior domains; first/last above-threshold crossings; censored if support reaches boundary.',
        'coarse_episode_windows':'Historical curated trough coordinates identify cycles/domains only. Raw maximum and side-minimum-to-peak boundaries recomputed; no historical peak coordinate is an outcome.',
        'envelope_vs_connected':'rise_q/fall_q are whole-episode envelope. connected_rise_q/connected_fall_q bracket selected maximum using same coarse side minima. Neither uses nearest pixel-scale local minima as pulse baseline.',
        'pairrule_caution':'Gene-specific windows are supplemented by identical EVE-defined windows for all genes. Common-window peaks at edges or without an identified homologous episode are ineligible; raw split-stripe peaks are anchor sensitivity, not sole primary biological reference.',
        'eve_centroid':'Integral of x*(y-local_min) divided by integral(y-local_min) over coarse cycle, exact for piecewise-linear raw profile.',
        'eve_relative_phase':'No extrapolation; fractional distance is descriptive and can change under nonlinear time-space maps. Ordinal before/after anchor retains meaning conditional on valid homologous anchors.',
        'eve_anchor_filters':'Stage4.1 cycle4 developing posterior pulse omitted from phase anchors. Split raw EVE peaks retained only as sensitivity; centroid and ALL prominent-lobe controls required. Kr4.3 halfheight is before centroid but not before every EVE2 lobe in every embryo, so this is not an anchor-independent order reversal.',
        'neural_exclusion':'Stage5.2 Kr and svb displayed for provenance but excluded from primary comparisons because likely neural expression was identified in manuscript (stronger for svb).',
        'extra_files':'EMB4 at4.1/4.3 separately flagged provenance_unverified; not included in n=3, SEM, historical cycle mapping, or primary inference.',
        'statistics':'Three biological embryos per stage. t95 intervals descriptive under normal-distribution assumption; exact two-sided sign test gives p=.25 for 3/3 concordance. Spatial pixels, landmarks, and thresholds are not replicates.'
    },indent=2))
    print('Spatial outputs:',len(ee),'episodes,',len(dd),'domain records,',len(phase_rows),'relative events.')
    print(aa[['stage','embryo','cycle','peak','centroid','multi_peak_flag','peak_anchor_primary_eligible']].to_string(index=False))

if __name__=='__main__':main()
