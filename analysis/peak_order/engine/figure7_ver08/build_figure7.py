"""Figure 7 ver08: submission gene palette; direct PDF and editable SVG.

The underlying curves, peak assignments and 20,000 biological bootstrap draws
are unchanged from version 05. Differences come from the stored draw-wise
contrasts, never by subtracting the bounds of two marginal peak intervals.
"""
import json
import numpy as np
import pandas as pd
from plot_core import *

F=VectorFigure(540,856)
MARKS=[]; CURVES=[]; DIFFS=[]
pr=NEW/'pairrule'; gp=NEW/'gap_peaks'
SP=pd.concat([pd.read_csv(pr/'spatial_peak_references.csv'),pd.read_csv(gp/'spatial_peak_references.csv')],ignore_index=True)
TP=pd.concat([pd.read_csv(pr/'temporal_peak_references.csv'),pd.read_csv(gp/'temporal_peak_references.csv')],ignore_index=True)
MASK=pd.read_csv(pr/'posterior_mask.csv')
PRC=pd.read_csv(pr/'all_temporal_pairwise_contrasts.csv')
GPC=pd.read_csv(gp/'temporal_pairwise_contrasts.csv')

def row_temporal(cohort,event):
    a=TP[(TP.cohort.str.lower()==cohort)&(TP.event_id==event)]
    assert len(a)==1, (cohort,event,len(a))
    return a.iloc[0]
def row_spatial(stage,emb,event):
    a=SP[(SP.stage==stage)&(SP.embryo==emb)&(SP.event_id==event)]
    assert len(a)==1, (stage,emb,event,len(a))
    return a.iloc[0]
def mask_lower(stage,emb):
    a=MASK[(MASK.stage==stage)&(MASK.embryo==emb)]
    for key in ['mask_ap_lo','verified_ap_lo']:
        if key in a and a[key].notna().any():return float(a[key].dropna().iloc[0])
    raise ValueError('No verified posterior region')
def valid_spatial(a):
    excluded=a.get('excluded_neural',False)
    return bool(a.spatial_resolved) and bool(a.in_posterior_mask) and not (pd.notna(excluded) and bool(excluded))

def contrast(panel,label,a,b,cohort=None,reverse=False,color='#333333'):
    if cohort=='pairrule':
        d=PRC[(PRC.event_a==a)&(PRC.event_b==b)]
        assert len(d)==1
        r=d.iloc[0]; est,lo,hi=map(float,[r.estimate,r.ci95_low,r.ci95_high])
        ca=cb='pairrule'; matched=True; src=pr/'all_temporal_pairwise_contrasts.csv'
    else:
        d=GPC[(GPC.event_a==a)&(GPC.event_b==b)]
        if cohort:d=d[(d.cohort_a.str.lower()==cohort)&(d.cohort_b.str.lower()==cohort)]
        assert len(d)==1, (a,b,cohort,len(d))
        r=d.iloc[0];est,lo,hi=map(float,[r.temporal_difference_b_minus_a,r.temporal_ci95_low,r.temporal_ci95_high])
        ca,cb=r.cohort_a,r.cohort_b;matched=bool(r.matched_temporal_channels);src=gp/'temporal_pairwise_contrasts.csv'
    if reverse:est,lo,hi=-est,-hi,-lo
    assert lo<=est<=hi
    out=dict(panel=panel,label=label,event_first=a if reverse else b,event_subtracted=b if reverse else a,
             cohort_first=ca if reverse else cb,cohort_subtracted=cb if reverse else ca,
             estimate=est,ci95_low=lo,ci95_high=hi,valid_bootstrap_fraction=float(r.valid_bootstrap_fraction),
             valid_bootstrap_n=int(r.valid_bootstrap_n),matched_temporal_channels=matched,
             units='ordinal_substage_steps',source=str(src),source_event_a=a,source_event_b=b,
             source_sign_reversed=reverse,color=color)
    DIFFS.append(out);return out

C={
 'A':[contrast('A','run - eve','eve_cycle3','run_cycle3','pairrule',color=COL['run']),
      contrast('A','odd - run*','run_cycle3','odd_cycle3','pairrule',color=COL['odd'])],
 'B':[contrast('B','svb - eve','svb_first','eve_cycle3','svb',True,COL['svb'])],
 'C':[contrast('C','gt1 - mlpt','mlpt_first','gt_peak1',color=COL['gt']),
      contrast('C','gt2 - mlpt','mlpt_first','gt_peak2',color=COL['gt'])],
 'D':[contrast('D','Kr - eve','kr_first','eve_cycle3','kr',True,COL['kr'])]
}

def marker(fx,fy,x,y,a,gene,domain,panel,stage='',embryo='',cohort=''):
    pos=float(a.estimate if domain=='time' else a.position)
    if not np.isfinite(pos):return
    if pos<x[0]-1e-9 or pos>x[-1]+1e-9:raise ValueError(f'Peak outside plotted curve: {a.event_id}')
    yy=float(np.interp(pos,x,y));split=domain=='space' and a.reference_kind=='split_domain_centre'
    lo=hi=np.nan
    if domain=='time':
        lo,hi=float(a.ci95_low),float(a.ci95_high)
        assert x[0]<=lo<=hi<=x[-1], (a.event_id,lo,hi,x[0],x[-1])
        F.line(fx(lo),fy(yy),fx(hi),fy(yy),COL[gene],.9)
        for xx in [lo,hi]:F.line(fx(xx),fy(yy)-2.5,fx(xx),fy(yy)+2.5,COL[gene],.9)
    if domain=='space' and float(a.lobe_hi)>float(a.lobe_lo):
        for xx in [float(a.lobe_lo),float(a.lobe_hi)]:
            F.marker(fx(xx),fy(float(np.interp(xx,x,y))),COL[gene],filled=False,r=2.1)
    F.marker(fx(pos),fy(yy),COL[gene],kind='diamond' if split else 'circle',filled=True)
    if gene=='gt':F.text(fx(pos)+4,fy(yy)-5,'1' if a.event_id=='gt_peak1' else '2',8,color=COL[gene])
    MARKS.append(dict(panel=panel,domain=domain,stage=stage,embryo=embryo,cohort=cohort,event_id=a.event_id,gene=gene,position=pos,display_y=yy,reference_kind=a.get('reference_kind','single_peak'),ci95_low=lo,ci95_high=hi,valid_bootstrap_fraction=float(a.valid_bootstrap_fraction) if domain=='time' else np.nan))

def key(x,y,genes):
    for gene in genes:
        F.line(x,y-2,x+11,y-2,COL[gene],1.1);F.text(x+15,y,'Kr' if gene=='kr' else gene,8.5,color=COL[gene],italic=True)
        x+=49 if gene!='mlpt' else 56

def difference_strip(letter,ox,oy):
    F.text(ox+26,oy+148,'Peak-time differences (95% intervals)',8,color='#333333')
    fx=lambda v:ox+103+(v+.5)/3.5*151
    F.line(fx(0),oy+154,fx(0),oy+194,'#999999',.6)
    rows=[164,182] if len(C[letter])==2 else [173]
    for d,yy in zip(C[letter],rows):
        y=oy+yy;lo,hi=d['ci95_low'],d['ci95_high']
        assert -.5<=lo<=hi<=3
        F.text(ox+26,y+2.5,d['label'],8,color='#333333')
        F.line(fx(lo),y,fx(hi),y,d['color'],1)
        for end in [lo,hi]:F.line(fx(end),y-2.5,fx(end),y+2.5,d['color'],.9)
        F.marker(fx(d['estimate']),y,d['color'],r=2.2)
    F.line(fx(-.5),oy+196,fx(3),oy+196,'#444444',.6)
    for v in [-.5,0,1,2,3]:
        F.line(fx(v),oy+196,fx(v),oy+199)
        F.text(fx(v),oy+209,str(v),8,anchor='middle')
    F.text(ox+141,oy+221,'Substage steps; positive = first named later',8,anchor='middle',color='#555555')

def draw_panel(letter,origin,title,outcome,series,time_xlim,space_xlim,cycles_only=False):
    ox,oy=origin;x=ox+42;w=214
    F.text(ox+9,oy+16,letter,12,bold=True)
    F.text(ox+26,oy+16,title,10,bold=True)
    key(ox+26,oy+32,[g for _,g,_ in series])
    F.text(ox+26,oy+46,'Time: '+outcome,8,color='#555555')
    ticks=[(i,LABELS[i]) for i in range(int(np.ceil(time_xlim[0])),int(np.floor(time_xlim[1]))+1)]
    if len(ticks)>7:ticks=ticks[::2]
    fx,fy=F.axes(x,oy+55,w,50,time_xlim,ticks,'Developmental stage',ylabel=True)
    for cohort,gene,events in series:
        xx,yy,src=temporal_curve(cohort,gene);lo,hi=time_xlim
        if gene=='eve' and not cycles_only:
            tr=row_temporal(cohort,events[0]);lo=max(lo,float(tr.window_lo));hi=min(hi,float(tr.window_hi))
        xc,yc=crop(xx,yy,lo,hi);F.curve(fx,fy,xc,yc,COL[gene])
        CURVES.append(dict(panel=letter,domain='time',gene=gene,cohort=cohort,source=src,plot_min=lo,plot_max=hi))
        for event in events:marker(fx,fy,xc,yc,row_temporal(cohort,event),gene,'time',letter,cohort=cohort)
    difference_strip(letter,ox,oy)
    F.text(ox+26,oy+234,'Space: stage 4.3, anterior to posterior',8,color='#555555')
    for emb,y0 in [(1,oy+247),(2,oy+296),(3,oy+345)]:
        ticks=[(v/100,str(v)) for v in [20,40,60,80,100] if space_xlim[0]*100<=v<=space_xlim[1]*100] if emb==3 else []
        if cycles_only and emb==3:ticks=[(.45,'45'),(.60,'60'),(.75,'75')]
        fx,fy=F.axes(x,y0,w,36,space_xlim,ticks,'Position (% embryo length)' if emb==3 else '',embryo=emb,ylabel=emb==2)
        for cohort,gene,events in series:
            xx,yy,src=spatial_curve('Stage4_3',emb,gene);lo=max(space_xlim[0],mask_lower('Stage4_3',emb));hi=space_xlim[1]
            if gene=='eve' or cycles_only:
                rr=row_spatial('Stage4_3',emb,events[0]);lo=max(lo,float(rr.family_window_lo));hi=min(hi,float(rr.family_window_hi))
            xc,yc=crop(xx,yy,lo,hi);F.curve(fx,fy,xc,yc,COL[gene])
            CURVES.append(dict(panel=letter,domain='space',gene=gene,stage='Stage4_3',embryo=emb,source=src,plot_min=lo,plot_max=hi))
            for event in events:
                rr=row_spatial('Stage4_3',emb,event)
                if valid_spatial(rr):marker(fx,fy,xc,yc,rr,gene,'space',letter,'Stage4_3',emb,cohort)

draw_panel('A',(0,0),'eve, run and odd: cycle 3','eve, then run, then odd',
           [('pairrule','eve',['eve_cycle3']),('pairrule','run',['run_cycle3']),('pairrule','odd',['odd_cycle3'])],(5,9.5),(.44,.78),True)
draw_panel('B',(270,0),'svb and eve cycle 3','eve, then svb',
           [('svb','svb',['svb_first']),('svb','eve',['eve_cycle3'])],(3,10),(.20,1.0))
draw_panel('C',(0,415),'gt and mlpt','separate cohorts, aligned by stage',
           [('gt','gt',['gt_peak1','gt_peak2']),('mlpt','mlpt',['mlpt_first'])],(4,12),(.20,1.0))
draw_panel('D',(270,415),'Kr and eve cycle 3','peak order unresolved',
           [('kr','kr',['kr_first']),('kr','eve',['eve_cycle3'])],(3,10),(.20,1.0))
F.text(270,835,'Circles: peaks. Open circles: eve daughters. Diamonds: split-domain centres.',8,anchor='middle',color='#555555')
F.text(270,848,'Horizontal bars: 95% timing intervals. *odd3 recovered in 87% of bootstrap draws.',8,anchor='middle',color='#555555')

pd.DataFrame(MARKS).to_csv(ROOT/'Fig7_peak_markers_ver08.csv',index=False)
pd.DataFrame(CURVES).to_csv(ROOT/'Fig7_curve_sources_ver08.csv',index=False)
pd.DataFrame(DIFFS).to_csv(ROOT/'Fig7_peak_time_differences_ver08.csv',index=False)
stats=F.save('Fig7_ver08')
print(json.dumps(stats,indent=2))
