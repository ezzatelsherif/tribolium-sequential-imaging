"""Publication figures from saved analyses; does not refit the model.

python make_manuscript_figures.py --numerical-only
python make_manuscript_figures.py --numerical-only --movie
Provide --original-figure PATH to reassemble the complete Figure 3 using an
original PDF-compatible Illustrator source, which is read without alteration.
"""
from pathlib import Path
import argparse, json, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import fitz

PACKAGE=Path(__file__).resolve().parent
AP=argparse.ArgumentParser()
AP.add_argument('--analysis-dir',type=Path,default=PACKAGE)
AP.add_argument('--original-figure',type=Path,default=None)
AP.add_argument('--output-dir',type=Path,default=PACKAGE/'figures')
AP.add_argument('--movie',action='store_true')
AP.add_argument('--numerical-only',action='store_true',help='Reproduce numerical panels, S3/S4 and optional movie without the original microscopy composite')
ARGS=AP.parse_args()
OUT=ARGS.output_dir.resolve();OUT.mkdir(parents=True,exist_ok=True)
ARGS.numerical_only=ARGS.numerical_only or ARGS.original_figure is None
ROOT=ARGS.analysis_dir.resolve(); sys.path.insert(0,str(ROOT))
from pulse_fit import PulseFit,load,GENES,X
from pulse_forward import pulse_values
from sweep import simulate,features,pulse

BLUE='#236A9E'; ORANGE='#BA5332'; COLORS=[BLUE,ORANGE]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
 'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.7,'savefig.facecolor':'white'})
PARTS=OUT/'_parts';PARTS.mkdir(exist_ok=True)
GEOM=json.loads((ROOT/'data/geometry.json').read_text())['genes']
PAIRS=load()
PDF_CACHE={}

def savefig(fig,name,final=False):
    dest=OUT if final else PARTS
    fig.savefig(dest/f'{name}.pdf')
    if final:
        fig.savefig(dest/f'{name}.svg')
        fig.savefig(dest/f'{name}.png',dpi=160)
    plt.close(fig)

def signal_data(gene):
    means=np.array([p.mean(0) for p in PAIRS[gene]])
    sem=np.array([p.std(0,ddof=1)/np.sqrt(len(p)) for p in PAIRS[gene]])
    scale=means.max(0)
    return means/scale,sem/scale

def fit_data(gene,shape='cosine',mode='mapped'):
    return json.loads((ROOT/'fitresults'/f'fit_{gene}_{shape}.json').read_text())[mode]

def standard_stage(ax,labels=True):
    ax.set_xlim(.6,24.4)
    ax.set_xticks(np.arange(1,25,3),[f'{j}.1' for j in range(1,9)] if labels else [])
    ax.tick_params(labelsize=8.4,length=3)
    if labels:ax.set_xlabel('eve-en stage',fontsize=9,labelpad=2)

def panel_b():
    t,i,e,_=simulate(GEOM['Kr'],.1,.1,.3,.5)
    i=i/i.max();e=e/e.max()
    pd.DataFrame({'time_over_D':t,'intron':i,'exon':e}).to_csv(OUT/'illustrative_profiles.csv',index=False)
    fig,axes=plt.subplots(2,1,figsize=(3.25,3.0))
    for y,c,lab in [(i,BLUE,'Intronic'),(e,ORANGE,'Exonic')]:
        axes[0].plot(t,y,color=c,lw=1.8,label=lab)
    axes[0].plot(t,pulse(t,.5),color='.65',ls=':',lw=1,label='Initiation')
    axes[0].set_xlim(0,1.8);axes[0].set_ylim(-.03,1.1)
    axes[0].set_xlabel('Time / D',labelpad=1)
    axes[0].legend(frameon=False,fontsize=8,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.3),handlelength=1.1,columnspacing=.8)
    x=np.linspace(-.5,1.3,800);snap=.5
    for y,c in [(i,BLUE),(e,ORANGE)]:axes[1].plot(x,np.interp(snap+x,t,y,left=0,right=0),color=c,lw=1.8)
    axes[1].axvline(0,color='#E19A20',lw=1.2)
    axes[1].annotate('',xy=(-.35,.75),xytext=(-.08,.75),arrowprops=dict(arrowstyle='->',color='.55',lw=1.5))
    axes[1].text(.03,.94,'Fixed position',color='#9C6200',fontsize=8,transform=axes[1].transAxes)
    axes[1].set_xlim(-.5,1.3);axes[1].set_ylim(-.03,1.1)
    axes[1].set_xlabel('Position / (vD)',labelpad=1)
    for ax in axes:
        ax.set_yticks([0,1]);ax.set_ylabel('Signal / peak',fontsize=9,labelpad=1)
        ax.tick_params(labelsize=8,length=3)
    fig.subplots_adjust(left=.15,right=.99,top=.87,bottom=.13,hspace=.72)
    savefig(fig,'panel_B')
    return t,i,e

def map_panel(gene):
    g=GEOM[gene]
    fig,ax=plt.subplots(figsize=(2.55,.8))
    end=g['transcript_end_bp']
    ax.plot([0,end],[.42,.42],color='.55',lw=.9)
    for a,b in g['exon_targets']:ax.add_patch(Rectangle((a,.29),b-a,.26,color=ORANGE,lw=0))
    for a,b in g['intron_targets']:ax.add_patch(Rectangle((a,.75),b-a,.12,color=BLUE,lw=0))
    ax.text(0,1.1,gene,fontsize=11,fontweight='bold',fontstyle='italic')
    ax.text(8950,1.05,f'{end/1000:.3f} kb',ha='right',fontsize=8,color='.35')
    ax.text(0,-.02,"5′",fontsize=8,va='top')
    ax.text(end,-.02,"3′",fontsize=8,va='top',ha='right')
    ax.plot([6500,8500],[-.23,-.23],color='.2',lw=1)
    ax.text(7500,-.27,'2 kb',fontsize=7.5,ha='center',va='top')
    ax.set_xlim(-80,9100);ax.set_ylim(-.65,1.45);ax.axis('off')
    fig.subplots_adjust(left=.025,right=.98,top=1,bottom=0)
    savefig(fig,f'map_{gene}')

def data_panel(gene):
    y,se=signal_data(gene);fit=fit_data(gene)
    dense=np.linspace(1,24,1151)
    pred=PulseFit(gene,dt=fit['dt']).predict(fit,dense)
    fig,ax=plt.subplots(figsize=(3.35,1.63))
    for c,color in enumerate(COLORS):
        ax.plot(dense,pred[:,c],color=color,lw=1.6)
        ax.errorbar(X,y[:,c],yerr=se[:,c],fmt='o',color=color,ms=2.9,lw=.65,capsize=1.3,zorder=3)
    ax.set_ylim(-.05,1.2);ax.set_yticks([0,.5,1])
    standard_stage(ax)
    if gene=='mlpt':ax.set_ylabel('Signal / channel maximum',fontsize=9)
    fig.subplots_adjust(left=.16,right=.99,top=.93,bottom=.28)
    savefig(fig,f'data_{gene}')

def offsets_panel():
    df=pd.read_csv(ROOT/'empirical/stage_lag_bootstrap.csv')
    order=[('KR',1),('MLPT',1),('MLPT',2),('SVB',1),('SVB',2)]
    fig,ax=plt.subplots(figsize=(3.75,2.62))
    ax.axvline(0,color='.5',ls='--',lw=.8)
    for j,(g,p) in enumerate(order):
        for phase,shift,color in [('rising',-.15,BLUE),('falling',.15,ORANGE)]:
            r=df[(df.gene==g)&(df.pulse==p)&(df.metric==phase)].iloc[0]
            ax.errorbar(r.estimate,j+shift,xerr=[[r.estimate-r.lo95],[r.hi95-r.estimate]],fmt='o',ms=4,lw=1.1,capsize=2,color=color)
    ax.set_yticks(range(5),['Kr','mlpt 1','mlpt 2','svb 1','svb 2'],fontstyle='italic',fontsize=9)
    ax.set_ylim(4.6,-.6);ax.set_xlim(-1.2,1.7)
    ax.set_xticks([-1,0,1]);ax.grid(axis='x',alpha=.12)
    ax.set_xlabel('Exon minus intron crossing\n(eve-en substage units)',fontsize=9)
    ax.legend(handles=[Line2D([],[],color=BLUE,marker='o',lw=1,label='Rising'),Line2D([],[],color=ORANGE,marker='o',lw=1,label='Falling')],ncol=2,frameon=False,fontsize=9,loc='upper center',bbox_to_anchor=(.5,1.2))
    fig.subplots_adjust(left=.2,right=.97,top=.85,bottom=.25)
    savefig(fig,'offsets')

def place(page,path,dest,clip=None):
    key=str(path)
    if key not in PDF_CACHE:PDF_CACHE[key]=fitz.open(path)
    src=PDF_CACHE[key]
    page.show_pdf_page(fitz.Rect(dest),src,0,clip=fitz.Rect(clip) if clip else None)

def text(page,x,y,s,size=10,bold=False,color=(0,0,0)):
    page.insert_text((x,y),s,fontsize=size,fontname='hebo' if bold else 'helv',color=color)

def composite_fig3():
    doc=fitz.open();page=doc.new_page(width=800,height=878)
    original=ARGS.original_figure
    text(page,15,20,'Figure 3',14,True)
    text(page,15,45,'A',14,True)
    place(page,original,(30,41,277,283),(180,392,675,876))
    text(page,30,299,'Selected stages and age anchors',9,color=(.3,.3,.3))
    text(page,298,45,'B',14,True)
    text(page,322,45,'Illustrative mapped-target model',10,True)
    place(page,PARTS/'panel_B.pdf',(294,51,568,285))
    text(page,322,299,'Kr; T/D = 0.1, S/D = 0.1, L/D = 0.3',9,color=(.3,.3,.3))
    text(page,601,45,'C',14,True)
    place(page,original,(615,62,782,148),(1290,493,1527,615))
    text(page,601,183,"C'",14,True)
    place(page,original,(615,204,782,277),(1290,741,1527,845))
    text(page,610,299,'Posterior measurement regions',9,color=(.3,.3,.3))
    text(page,15,328,'D',14,True)
    text(page,40,328,'Mapped targets and representative embryos',10,True)
    text(page,267,348,'Posterior temporal profiles',10,True)
    clips={'Kr':(234,1007,635,1103),'mlpt':(234,1205,635,1302),'svb':(234,1404,635,1504)}
    for row,gene in enumerate(GENES):
        top=368+row*154
        place(page,PARTS/f'map_{gene}.pdf',(32,top,233,top+68))
        place(page,original,(32,top+77,233,top+128),clips[gene])
        place(page,PARTS/f'data_{gene}.pdf',(245,top-8,505,top+126))
    text(page,42,349,'Intronic target',8,color=(.137,.416,.620))
    text(page,133,349,'Exonic target',8,color=(.729,.325,.196))
    text(page,267,364,'Points: mean +/- SEM; lines: primary fit',8.5,color=(.3,.3,.3))
    text(page,32,838,'Targets share a common scale; zero is the aligned 5-prime target boundary.',8.5,color=(.3,.3,.3))
    text(page,267,815,'Intronic',9,color=(.137,.416,.620))
    text(page,350,815,'Exonic',9,color=(.729,.325,.196))
    text(page,267,854,'Equal spacing of eve-en substages',8.5,color=(.3,.3,.3))
    text(page,515,328,'E',14,True)
    place(page,original,(516,342,788,643),(989,925,1580,1570))
    text(page,515,672,'F',14,True)
    text(page,539,672,'Observed half-height offsets',10,True)
    place(page,PARTS/'offsets.pdf',(513,690,793,858))
    text(page,539,872,'Bars: nominal 95% paired-bootstrap intervals',8.3,color=(.3,.3,.3))
    doc.set_metadata({'title':'Figure 3. Intronic and exonic signal profiles','subject':'Revised analysis; original microscopy preserved','creator':'make_figures.py'})
    doc.save(OUT/'Figure_3.pdf',garbage=3,deflate=True)
    page.get_pixmap(matrix=fitz.Matrix(1.6,1.6),alpha=False).save(OUT/'Figure_3.png')
    (OUT/'Figure_3.svg').write_text(page.get_svg_image(text_as_path=False))
    doc.close()

def supplementary_s3():
    df=pd.read_csv(ROOT/'sweep_results/sweep_results.csv')
    grid=(.01,.03,.1,.3,1.,3.)
    fig=plt.figure(figsize=(10.2,10.4))
    gs=fig.add_gridspec(3,3,left=.08,right=.89,top=.90,bottom=.115,hspace=.50,wspace=.25,height_ratios=[1,1,1.15])
    for col,gene in enumerate(GENES):
        d=df[df.gene==gene]
        for row,ph in enumerate(('rise','fall')):
            ax=fig.add_subplot(gs[row,col])
            a=np.array([[np.mean(d[(d['T']==T)&(d.decay_mean==L)][f'{ph}_50_offset']>.01) for L in grid] for T in grid])
            im=ax.imshow(a,origin='lower',vmin=0,vmax=1,cmap='Blues')
            for j in range(6):
                for k in range(6):ax.text(k,j,str(round(a[j,k]*18)),ha='center',va='center',fontsize=8,color='white' if a[j,k]>.57 else 'black')
            ax.set_xticks(range(6),[str(x) for x in grid],fontsize=8)
            ax.set_yticks(range(6),[str(x) for x in grid],fontsize=8)
            ax.set_title(f'{gene}: {"rising" if ph=="rise" else "falling"}',fontstyle='italic',fontsize=11)
            if col==0:ax.set_ylabel('Traversal duration T / D')
            if row==1:ax.set_xlabel('Mean decay wait L / D')
        ax=fig.add_subplot(gs[2,col])
        for r,c,lab in [(.25,'#0072B2','Early peak'),(.5,'#D55E00','Symmetric'),(.75,'#009E73','Late peak')]:
            q=d[d.peak_fraction==r]
            ax.scatter(q.rise_50_offset,q.fall_50_offset,s=10,color=c,alpha=.5,linewidths=0,label=lab)
        ax.axvspan(-.01,.01,color='.65',alpha=.3);ax.axhspan(-.01,.01,color='.65',alpha=.3)
        ax.axvline(0,color='.45',lw=.7);ax.axhline(0,color='.45',lw=.7)
        ax.set_xlim(-1.05,1.55);ax.set_ylim(-.13,3.75)
        ax.set_title(gene,fontstyle='italic',fontsize=11)
        ax.set_xlabel('Rising half-height offset / D')
        if col==0:ax.set_ylabel('Falling half-height offset / D')
        if col==2:ax.legend(fontsize=8,frameon=False,loc='upper left')
    cb=fig.colorbar(im,cax=fig.add_axes([.925,.45,.018,.41]));cb.set_label('Fraction with intronic lead > 0.01 D',fontsize=9)
    fig.text(.025,.925,'A',fontsize=15,fontweight='bold')
    fig.text(.025,.36,'B',fontsize=15,fontweight='bold')
    fig.suptitle('S3 Fig. Dependence of signal ordering on biochemical parameters',y=.985,fontsize=13,fontweight='bold')
    fig.text(.5,.948,'1,944 hypothetical scenarios; 648 per gene. Each map cell contains 18 scenarios (6 processing waits x 3 pulse shapes).',ha='center',fontsize=9)
    fig.text(.5,.019,'Positive offsets: the intronic profile crosses first. D is the hypothetical initiation-pulse duration.\nCounts summarize the specified grid and are not biological probabilities. Grey bands mark offsets within +/-0.01 D.',ha='center',fontsize=9,linespacing=1.4)
    savefig(fig,'Figure_S3',final=True)

def supplementary_s4():
    fig,axes=plt.subplots(3,3,figsize=(10.4,8.7),gridspec_kw={'height_ratios':[2.1,1.2,1.]})
    dense=np.linspace(1,24,1151)
    for col,gene in enumerate(GENES):
        y,se=signal_data(gene)
        for shape,style in [('cosine','-'),('gaussian','--')]:
            fit=fit_data(gene,shape)
            pred=PulseFit(gene,shape=shape,dt=fit['dt']).predict(fit,dense)
            for c,color in enumerate(COLORS):axes[0,col].plot(dense,pred[:,c],style,color=color,lw=1.6)
            init=sum(a*pulse_values(dense,q['center'],q['rise'],q['fall'],shape) for a,q in zip(fit['amplitudes'],fit['parameters']['pulses']))
            axes[2,col].plot(dense,init/init.max(),style,color='.25',lw=1.5)
        primary=fit_data(gene)
        for c,color in enumerate(COLORS):
            axes[0,col].errorbar(X,y[:,c],yerr=se[:,c],fmt='o',color=color,ms=3,lw=.7,capsize=1.4)
            axes[1,col].errorbar(X,y[:,c]-np.array(primary['pred'])[:,c],yerr=se[:,c],fmt='o-',color=color,ms=3,lw=.7,capsize=1.4)
        axes[0,col].set_title(gene,fontsize=12,fontstyle='italic',fontweight='bold',loc='left')
        axes[0,col].set_ylim(-.08,1.23);axes[1,col].axhline(0,color='.5',ls=':',lw=1)
        axes[2,col].set_ylim(-.03,1.13)
        for row in range(3):standard_stage(axes[row,col],labels=(row==2));axes[row,col].grid(axis='y',alpha=.1)
    axes[0,0].set_ylabel('Signal / observed channel maximum')
    axes[1,0].set_ylabel('Observed minus\nprimary fit',fontsize=9)
    axes[2,0].set_ylabel('Initiation /\nown peak',fontsize=9)
    fig.legend(handles=[Line2D([],[],color=BLUE,marker='o',lw=0,label='Intronic'),Line2D([],[],color=ORANGE,marker='o',lw=0,label='Exonic'),Line2D([],[],color='.25',lw=1.5,label='Primary: raised cosine'),Line2D([],[],color='.25',ls='--',lw=1.5,label='Sensitivity: Gaussian')],ncol=4,frameon=False,bbox_to_anchor=(.5,.96),loc='upper center',fontsize=9)
    for y,lab in [(.883,'A'),(.45,'B'),(.265,'C')]:fig.text(.006,y,lab,fontsize=14,fontweight='bold')
    fig.suptitle('S4 Fig. Constrained fits, residuals, and shared transcription inputs',y=.99,fontsize=13,fontweight='bold')
    fig.text(.5,.018,'Points: stage means +/- 1 SEM (4-5 paired specimens). Both channels fitted jointly with shared kinetics across pulses.\nSubstages are equally spaced. Inferred durations are in stage units, not measured biochemical lifetimes.',ha='center',fontsize=9,linespacing=1.4)
    fig.tight_layout(rect=(.02,.065,1,.9),h_pad=1.5,w_pad=1.6)
    savefig(fig,'Figure_S4',final=True)

def movie(t,i,e):
    from matplotlib.animation import FuncAnimation,FFMpegWriter
    plt.rcParams.update({'font.size':11})
    fig,axes=plt.subplots(1,2,figsize=(9,4.6))
    fig.suptitle('Mapped-target model: temporal persistence and a travelling profile',fontsize=13)
    axes[0].plot(t,i,color=BLUE,lw=2,label='Intronic');axes[0].plot(t,e,color=ORANGE,lw=2,label='Exonic')
    axes[0].plot(t,pulse(t,.5),color='.6',ls=':',lw=1.2,label='Initiation')
    axes[0].set_xlim(0,1.8);axes[0].set_xlabel('Time / D');axes[0].set_ylabel('Signal / channel peak')
    axes[0].legend(frameon=False,fontsize=9)
    cursor=axes[0].axvline(0,color='#D8941A',lw=1.2)
    x=np.linspace(-1.8,1.8,900)
    li,=axes[1].plot([],[],color=BLUE,lw=2);le,=axes[1].plot([],[],color=ORANGE,lw=2)
    axes[1].set_xlim(-1.8,1.8);axes[1].set_xlabel('Position / (vD)')
    axes[1].axvline(0,color='#D8941A',lw=1.2)
    axes[1].annotate('Propagation',xy=(-1.4,.83),xytext=(-.65,.83),arrowprops=dict(arrowstyle='->',color='.4'),fontsize=9,ha='center',va='center')
    axes[1].set_title('Same response at shifted local times',fontsize=10)
    time_label=fig.text(.5,.12,'',ha='center',fontsize=11)
    for ax in axes:ax.set_ylim(-.03,1.1)
    fig.text(.5,.025,'Illustrative Kr geometry: T/D = 0.1, S/D = 0.1, L/D = 0.3. Fixed normalization throughout.\nD is hypothetical pulse duration; v is an arbitrary constant propagation speed.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.075,right=.985,top=.84,bottom=.27,wspace=.23)
    def update(z):
        cursor.set_xdata([z,z]);li.set_data(x,np.interp(z+x,t,i,left=0,right=0));le.set_data(x,np.interp(z+x,t,e,left=0,right=0))
        time_label.set_text(f'Time / D = {z:.2f}')
        return cursor,li,le,time_label
    anim=FuncAnimation(fig,update,frames=np.linspace(0,1.8,144),interval=50,blit=False)
    anim.save(OUT/'Movie_S1.mp4',writer=FFMpegWriter(fps=20,codec='libx264',extra_args=['-pix_fmt','yuv420p']),dpi=120)
    update(.5);fig.savefig(OUT/'Movie_S1_preview.png',dpi=140)
    plt.close(fig)

if __name__=='__main__':
    t,i,e=panel_b();offsets_panel()
    for gene in GENES:map_panel(gene);data_panel(gene)
    if not ARGS.numerical_only:composite_fig3()
    supplementary_s3();supplementary_s4()
    if ARGS.movie:movie(t,i,e)
    print('Created numerical panels, Figure_S3, Figure_S4, and requested composite/movie in',OUT)
