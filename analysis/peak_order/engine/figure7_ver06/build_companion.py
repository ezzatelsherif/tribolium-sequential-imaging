"""All-stage companion with peak-time intervals and direct temporal contrasts."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from plot_core import *

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
pr=NEW/'pairrule';gp=NEW/'gap_peaks'
TP=pd.concat([pd.read_csv(pr/'temporal_peak_references.csv'),pd.read_csv(gp/'temporal_peak_references.csv')],ignore_index=True)
SP=pd.concat([pd.read_csv(pr/'spatial_peak_references.csv'),pd.read_csv(gp/'spatial_peak_references.csv')],ignore_index=True)
M=pd.read_csv(pr/'posterior_mask.csv')
PAGES=[];MARKS=[];BARS=[]
pdf=PdfPages(ROOT/'Fig7_Peak_Order_Companion_ver06.pdf')

def finish(fig,name):
    pdf.savefig(fig);fig.savefig(ROOT/(name+'.png'),dpi=140);plt.close(fig);PAGES.append(name)
def note(fig,s):fig.text(.08,.025,s,fontsize=8,color='#555555',va='bottom')
def mark(ax,x,y,row,gene,domain,page,stage='',embryo='',cohort=''):
    pos=float(row.estimate if domain=='time' else row.position)
    if not np.isfinite(pos) or not x[0]<=pos<=x[-1]:return
    yy=float(np.interp(pos,x,y));split=domain=='space' and row.reference_kind=='split_domain_centre'
    if split:
        for xx in [row.lobe_lo,row.lobe_hi]:
            ax.plot(xx*100,np.interp(xx,x,y),'o',ms=4,mfc='white',mec=COL[gene],mew=.8,zorder=4)
    if domain=='time' and pd.notna(row.ci95_low) and pd.notna(row.ci95_high):
        lo,hi=float(row.ci95_low),float(row.ci95_high)
        ax.hlines(yy,lo,hi,color=COL[gene],lw=1.15,zorder=4)
        ax.plot([lo,hi],[yy,yy],'|',color=COL[gene],ms=6,mew=1.0,zorder=4)
        BARS.append(dict(page=page,panel=ax.get_gid(),cohort=cohort,gene=gene,event_id=row.event_id,estimate=pos,ci95_low=lo,ci95_high=hi,display_y=yy,valid_bootstrap_fraction=row.valid_bootstrap_fraction,interval_type='marginal_peak_position_percentile95',source=str(pr/'temporal_peak_references.csv' if cohort=='pairrule' else gp/'temporal_peak_references.csv')))
    ax.plot(pos if domain=='time' else pos*100,yy,'D' if split else 'o',ms=4,mfc=COL[gene],mec='white',mew=.45,zorder=5)
    MARKS.append(dict(page=page,domain=domain,stage=stage,embryo=embryo,cohort=cohort,gene=gene,event_id=row.event_id,position=pos,display_y=yy,reference_kind=row.get('reference_kind','peak')))
    if gene=='eve':
        cycle=int(row.cycle);ax.annotate(str(cycle),(pos if domain=='time' else pos*100,yy),xytext=(0,6),textcoords='offset points',ha='center',fontsize=7.5,color=COL[gene])
    elif gene=='gt':
        ax.annotate('1' if row.event_id=='gt_peak1' else '2',(pos if domain=='time' else pos*100,yy),xytext=(5,4),textcoords='offset points',fontsize=7.5,color=COL[gene])
def temporal_panel(ax,cohort,genes,title,page):
    ax.set_gid(title)
    for gene in genes:
        xx,yy,_=temporal_curve(cohort,gene);x,y=crop(xx,yy,0,15)
        ax.plot(x,y,lw=1.2,color=COL[gene],label='Kr' if gene=='kr' else gene)
        rows=TP[(TP.cohort.str.lower()==cohort)&(TP.gene.str.lower()==gene)]
        rows=rows[(rows.cycle.isna())|(rows.cycle.between(2,5))]
        for _,row in rows.iterrows():mark(ax,x,y,row,gene,'time',page,cohort=cohort)
    ax.set_xlim(0,15);ax.set_ylim(-.04,1.18);ax.set_yticks([0,.5,1]);ax.set_ylabel('Relative expression')
    ax.set_xticks(range(0,16),LABELS[:16],fontsize=8)
    ax.set_title(title,loc='left',weight='bold',fontsize=10);ax.legend(frameon=False,ncols=len(genes),loc='upper left',fontsize=8)
    ax.set_xlabel('Developmental stage (ordered substages)')

fig,axs=plt.subplots(3,1,figsize=(8.3,11.7));fig.subplots_adjust(left=.10,right=.97,top=.93,bottom=.10,hspace=.53)
fig.suptitle('Temporal expression: matched measurement cohorts',x=.10,ha='left',fontsize=14,weight='bold')
for ax,args in zip(axs,[('pairrule',['eve','run','odd'],'eve, run and odd measured together'),('kr',['kr','eve'],'Kr and its matched eve reference'),('svb',['svb','eve'],'svb and its matched eve reference')]):temporal_panel(ax,*args,'Companion_01_temporal_matched')
note(fig,'Filled circles: smoothed peak estimates. Horizontal bars: 95% biological-bootstrap intervals for peak times.\nNumbers identify temporal eve cycles. Intervals are conditional on stages, smoothing and valid peak recovery.\nProfiles are scaled separately. These fixed-embryo reconstructions are not longitudinal trajectories.')
finish(fig,'Companion_01_temporal_matched')

fig,axs=plt.subplots(3,1,figsize=(8.3,11.7));fig.subplots_adjust(left=.10,right=.97,top=.93,bottom=.10,hspace=.53)
fig.suptitle('Temporal expression: gt and mlpt cohorts',x=.10,ha='left',fontsize=14,weight='bold')
temporal_panel(axs[0],'gt',['gt','eve'],'gt and its matched eve reference','Companion_02_temporal_gt_mlpt')
temporal_panel(axs[1],'mlpt',['mlpt','eve'],'mlpt and its matched eve reference','Companion_02_temporal_gt_mlpt')
ax=axs[2]
ax.set_gid('gt and mlpt aligned by assigned stage: separate cohorts')
for cohort,gene in [('gt','gt'),('mlpt','mlpt')]:
    xx,yy,_=temporal_curve(cohort,gene);x,y=crop(xx,yy,0,15);ax.plot(x,y,lw=1.2,color=COL[gene],label=gene)
    for _,row in TP[(TP.cohort.str.lower()==cohort)&(TP.gene==gene)].iterrows():mark(ax,x,y,row,gene,'time','Companion_02_temporal_gt_mlpt',cohort=cohort)
ax.set_xlim(0,15);ax.set_ylim(-.04,1.18);ax.set_xticks(range(0,16),LABELS[:16],fontsize=8);ax.set_ylabel('Relative expression');ax.set_xlabel('Developmental stage (ordered substages)')
ax.set_yticks([0,.5,1])
ax.set_title('gt and mlpt aligned by assigned stage: separate cohorts',loc='left',weight='bold',fontsize=10);ax.legend(frameon=False,ncols=2,loc='upper left',fontsize=8)
note(fig,'The bottom panel aligns independent cohorts by assigned stage; it does not create matched embryos.\nHorizontal bars show 95% bootstrap intervals for peak times, conditional on staging, smoothing and recovery.\nNumbers 1/2 identify gt peaks; eve numbers identify cycles. Later expression remains visible for context.')
finish(fig,'Companion_02_temporal_gt_mlpt')

groups=[(['eve','run','odd'],'eve / run / odd'),(['kr','eve'],'Kr / eve'),(['svb','eve'],'svb / eve'),(['gt','mlpt','eve'],'gt / mlpt / eve')]
for si,stage in enumerate(STAGES):
    page=f'Companion_{si+3:02d}_{stage}'
    fig,axs=plt.subplots(4,3,figsize=(8.3,11.7),sharex=True);fig.subplots_adjust(left=.085,right=.98,top=.86,bottom=.12,hspace=.63,wspace=.19)
    stlabel=stage[5:].replace('_','.')
    fig.suptitle(f'Spatial expression at stage {stlabel}',x=.085,ha='left',fontsize=14,weight='bold')
    complete=bool(M[M.stage==stage].mask_complete.all())
    fig.text(.085,.925,'Three most posterior eve domain families' if complete else 'Verified two posterior eve families; third-family boundary unresolved',fontsize=9,color='#555555')
    fig.legend([Line2D([0],[0],color=COL[g],lw=1.3) for g in ['eve','run','odd','kr','svb','gt','mlpt']],['eve','run','odd','Kr','svb','gt','mlpt'],loc='upper left',bbox_to_anchor=(.085,.912),ncols=7,frameon=False,fontsize=8)
    mm=M[M.stage==stage]
    lo=float(mm.verified_ap_lo.min())
    for row,(genes,title) in enumerate(groups):
        for j,emb in enumerate([1,2,3]):
            ax=axs[row,j];m=mm[mm.embryo==emb];cut=float(m.verified_ap_lo.iloc[0])
            for gene in genes:
                x,y,_=spatial_curve(stage,emb,gene);xc,yc=crop(x,y,cut,1.)
                if gene=='svb' and stage=='Stage5_2':
                    ax.set_title('svb excluded: neural expression',loc='left',fontsize=8,color=COL[gene]);continue
                ax.plot(xc*100,yc,lw=1.0,color=COL[gene],label=gene)
                rr=SP[(SP.stage==stage)&(SP.embryo==emb)&(SP.gene==gene)&SP.spatial_resolved&SP.in_posterior_mask]
                for _,a in rr.iterrows():
                    excluded=a.get('excluded_neural',False)
                    if pd.notna(excluded) and bool(excluded):continue
                    mark(ax,xc,yc,a,gene,'space',page,stage,emb)
                if gene=='kr' and stage=='Stage5_2':ax.set_title('Kr peak is anterior to this region',loc='left',fontsize=8,color=COL[gene])
            ax.set_xlim(lo*100,100);ax.set_ylim(-.04,1.20);ax.set_yticks([0,1]);ax.tick_params(labelsize=8)
            if row==0:ax.set_title(f'Embryo {emb}',fontsize=10,weight='bold')
            if j==0:ax.set_ylabel(title+'\nRelative expression',fontsize=8)
            if row==3:ax.set_xlabel('Anterior to posterior\n(% embryo length)',fontsize=8)
    note(fig,'Filled circles: peaks. Open circles: two eve daughter peaks. Diamonds: their midpoint (domain centre).\nEve numbers identify temporal cycles, not posterior rank. Each embryo has its own posterior cutoff.\nCurves are smoothed before cropping and scaled separately over the full profile; embryos remain separate.')
    finish(fig,page)

# Direct contrasts preserve covariance in matched channels; they are not differences
# between separately reported marginal interval endpoints.
PRC=pd.read_csv(pr/'all_temporal_pairwise_contrasts.csv')
GPC=pd.read_csv(gp/'temporal_pairwise_contrasts.csv')
CONTRASTS=[]
def contrast(plus,minus,label,group):
    data=PRC if plus.split('_')[0] in ['eve','run','odd'] and minus.split('_')[0] in ['eve','run','odd'] else GPC
    rr=data[((data.event_a==plus)&(data.event_b==minus))|((data.event_a==minus)&(data.event_b==plus))]
    assert len(rr)==1,(plus,minus,len(rr))
    r=rr.iloc[0]; sign=1 if r.event_b==plus else -1
    est=float(r.get('estimate',r.get('temporal_difference_b_minus_a')))
    lo=float(r.get('ci95_low',r.get('temporal_ci95_low')));hi=float(r.get('ci95_high',r.get('temporal_ci95_high')))
    lo,hi=(lo,hi) if sign==1 else (-hi,-lo)
    CONTRASTS.append(dict(plus_event=plus,minus_event=minus,label=label,group=group,estimate=sign*est,ci95_low=lo,ci95_high=hi,valid_bootstrap_fraction=float(r.valid_bootstrap_fraction),valid_bootstrap_n=int(r.valid_bootstrap_n),interval_type='direct_difference_percentile95_conditional_on_valid_features',source=str(pr/'all_temporal_pairwise_contrasts.csv' if data is PRC else gp/'temporal_pairwise_contrasts.csv')))

for a,b,l,g in [
    ('run_cycle3','eve_cycle3','run3 - eve3','Matched temporal channels'),
    ('odd_cycle3','run_cycle3','odd3 - run3 *','Matched temporal channels'),
    ('svb_first','eve_cycle3','svb - eve3','Matched temporal channels'),
    ('kr_first','eve_cycle3','Kr - eve3','Matched temporal channels'),
    ('odd_cycle2','run_cycle2','odd2 - run2','Matched temporal channels'),
    ('gt_peak1','mlpt_first','gt1 - mlpt','Independent temporal cohorts'),
    ('gt_peak2','mlpt_first','gt2 - mlpt','Independent temporal cohorts'),
    ('gt_peak1','run_cycle3','gt1 - run3','Independent temporal cohorts')]:contrast(a,b,l,g)
DC=pd.DataFrame(CONTRASTS)
fig=plt.figure(figsize=(8.3,11.7))
fig.suptitle('How far apart are the temporal peaks?',x=.08,y=.96,ha='left',fontsize=14,weight='bold')
fig.text(.08,.91,'Each point is the first named peak time minus the second.\nPositive values mean that the first peak occurs later.\nBars show direct 95% bootstrap intervals for that difference.\nNumbers denote expression cycles (eve/run/odd) or assigned peaks (gt).',fontsize=10,va='top',linespacing=1.5)
for group,box in [('Matched temporal channels',[.29,.48,.43,.30]),('Independent temporal cohorts',[.29,.19,.43,.20])]:
    ax=fig.add_axes(box);d=DC[DC.group==group].reset_index(drop=True)
    ax.axvline(0,color='#999999',lw=.8,ls='--')
    for i,row in d.iterrows():
        col='#1767a5' if group.startswith('Matched') else '#a56d16'
        ax.hlines(i,row.ci95_low,row.ci95_high,color=col,lw=1.5)
        ax.plot([row.ci95_low,row.ci95_high],[i,i],'|',color=col,ms=7,mew=1.1)
        ax.plot(row.estimate,i,'o',ms=5,color=col)
        ax.text(1.04,i,f'{row.estimate:+.3f}\n[{row.ci95_low:+.3f}, {row.ci95_high:+.3f}]',transform=ax.get_yaxis_transform(),va='center',fontsize=8,color='#444444')
    ax.set_title(group,loc='left',fontsize=11,weight='bold',pad=14)
    ax.set_yticks(range(len(d)),d.label,fontsize=10);ax.set_ylim(len(d)-.5,-.5);ax.set_xlim(-.5,3.0);ax.set_xticks([-.5,0,.5,1,1.5,2,2.5,3]);ax.tick_params(axis='y',length=0)
    ax.set_xlabel('Peak-time difference (ordered-substage intervals)',fontsize=9)
    ax.spines['left'].set_visible(False)
fig.text(.08,.105,'* The odd3 peak was recovered in 86.99% of bootstrap draws. Its interval is conditional on recovery.\nAll other displayed contrasts have 100% valid draws. The intervals quantify biological sampling under\nthe chosen stages and smoother; they do not include unknown age or cross-cohort staging offsets.\nThese are pointwise, unadjusted intervals. Spatial embryos remain separate and have no pixel-based error bars.',fontsize=8,color='#555555',va='top',linespacing=1.5)
finish(fig,'Companion_07_temporal_peak_differences')
pdf.close()
(ROOT/'companion_pages.json').write_text(json.dumps(PAGES,indent=2))
pd.DataFrame(MARKS).to_csv(ROOT/'Companion_peak_markers_ver06.csv',index=False)
pd.DataFrame(BARS).to_csv(ROOT/'Companion_temporal_peak_intervals_ver06.csv',index=False)
DC.to_csv(ROOT/'Companion_temporal_peak_differences_ver06.csv',index=False)
print('Created',len(PAGES),'pages')
