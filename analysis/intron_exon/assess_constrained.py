"""Focused conditional lack-of-fit assessment and traversal profiles.

This is a model-centered paired residual bootstrap, not predictive validation
and not an estimate of parameter identifiability. Channel normalizers stay fixed.
Run profiles before finalized fits, then null after central fits are finalized.
"""
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from pulse_fit import PulseFit, GENES, load, serialize

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'assessment'
TRAVERSALS=[.005,.03,.1,.3,1.,3.,6.]

def source(gene,dt=None):
    path=ROOT/'fitresults'/f'fit_{gene}_cosine.json'
    fits=json.loads(path.read_text())
    pairs=load()[gene]
    means=np.array([p.mean(0) for p in pairs]);scale=means.max(0)
    y=means/scale
    selected_dt=dt if dt is not None else fits['mapped']['dt']
    fitter=PulseFit(gene,'mapped','cosine',dt=selected_dt)
    return fitter,fits,pairs,scale,y,path

def fixed_traversal(fitter,y,T,starts,max_nfev=220):
    lo,hi=np.array(fitter.bounds).T
    j=fitter.names.index('T');keep=np.arange(len(lo))!=j
    fixed=np.log(T);candidates=[];diagnostics=[]
    def expand(z):
        theta=np.empty(len(lo));theta[keep]=z;theta[j]=fixed
        return theta
    for s in starts:
        theta=fitter.pack(s);theta[j]=fixed
        init=fitter.solve(theta,y);candidates.append(init)
        opt=least_squares(lambda z:fitter.solve(expand(z),y)['residual'].ravel(),
            theta[keep],bounds=(lo[keep],hi[keep]),max_nfev=max_nfev,
            ftol=1e-10,xtol=1e-9,gtol=1e-8,diff_step=1e-4)
        f=fitter.solve(expand(opt.x),y);candidates.append(f)
        diagnostics.append(dict(cost=f['cost'],nfev=opt.nfev,
            success=bool(opt.success),optimality=float(opt.optimality),
            message=opt.message))
    best=min(candidates,key=lambda f:f['cost']);best['optimization']=diagnostics
    return best

def profile(gene,dt=None):
    fitter,fits,pairs,scale,y,path=source(gene,dt)
    p=fits['mapped']['parameters'];c=dict(fits['collapsed']['parameters'],T=.005)
    starts=[p,c]
    alt=json.loads(json.dumps(p));alt['splice']=.05;alt['decay']=1.
    starts.append(alt)
    alt=json.loads(json.dumps(p));alt['splice']=1.;alt['decay']=.05
    starts.append(alt)
    rows=[];results=[];t0=time.time()
    for T in TRAVERSALS:
        f=fixed_traversal(fitter,y,T,starts)
        starts=starts[:4]+[f['parameters']]
        results.append(f)
        r=dict(gene=gene,T=T,SSE=f['cost'],dt=fitter.dt,
            successful_starts=sum(d['success'] for d in f['optimization']),
            total_starts=len(f['optimization']))
        r.update({k:v for k,v in f['parameters'].items() if k!='pulses'})
        for j,q in enumerate(f['parameters']['pulses']):
            r.update({f'{k}{j+1}':v for k,v in q.items()})
        rows.append(r)
        print('profile',gene,T,f['cost'],'seconds',round(time.time()-t0,1),flush=True)
    # A descending continuation helps detect path dependence between local basins.
    descending=results[-1]['parameters']
    for j in range(len(results)-2,-1,-1):
        f=fixed_traversal(fitter,y,TRAVERSALS[j],[descending,results[j]['parameters']])
        if f['cost']<results[j]['cost']:
            results[j]=f;rows[j].update(SSE=f['cost'])
            rows[j].update(successful_starts=sum(d['success'] for d in f['optimization']),
                total_starts=len(f['optimization']))
            rows[j].update({k:v for k,v in f['parameters'].items() if k!='pulses'})
            for pi,q in enumerate(f['parameters']['pulses']):rows[j].update({f'{k}{pi+1}':v for k,v in q.items()})
        descending=f['parameters']
    pd.DataFrame(rows).to_csv(OUT/f'traversal_profile_{gene}.csv',index=False)
    (OUT/f'traversal_profile_{gene}.json').write_text(json.dumps(serialize(dict(
        gene=gene,dt=fitter.dt,source_file=str(path.relative_to(ROOT)),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        fits=results,elapsed_seconds=time.time()-t0)),indent=2))
    best=min(results,key=lambda f:f['cost'])
    print('PROFILE COMPLETE',gene,'best',best['cost'],'central',fits['mapped']['cost'],flush=True)

def null(gene,B=199,dt=None):
    fitter,fits,pairs,scale,y,path=source(gene,dt)
    central=fitter.solve(fitter.pack(fits['mapped']['parameters']),y)
    if abs(central['cost']-fits['mapped']['cost'])>1e-7:
        raise ValueError('Stored fit cost does not match current dt and model; finalize central fit first.')
    # Choose a second plausible biochemical basin from the completed fixed-T profile.
    profile_path=OUT/f'traversal_profile_{gene}.json'
    prof=json.loads(profile_path.read_text())['fits']
    p0=central['parameters'];theta0=fitter.pack(p0)
    ranked=sorted(prof,key=lambda f:f['cost'])
    alternatives=[f for f in ranked if abs(np.log(f['parameters']['T']/p0['T']))>np.log(3)
        and f['parameters']['T']>=.3]
    alt=(alternatives or ranked)[0]['parameters']
    starts=[p0,alt]
    rng=np.random.default_rng(2026091300+GENES.index(gene))
    residuals=[(p-p.mean(0))/scale*np.sqrt(len(p)/(len(p)-1)) for p in pairs]
    rows=[];allfits=[];t0=time.time()
    for b in range(B):
        pseudo=central['pred']+np.array([r[rng.integers(0,len(r),size=len(r))].mean(0) for r in residuals])
        f=fitter.fit(pseudo,starts=starts,global_search=False)
        diagnostics=f['optimization']
        retried=False
        # A limit-hit in every start is a numerical failure requiring a longer refit.
        if not any(d['success'] for d in diagnostics):
            lo,hi=np.array(fitter.bounds).T
            opt=least_squares(lambda z:fitter.solve(z,pseudo)['residual'].ravel(),
                f['theta'],bounds=(lo,hi),max_nfev=500,ftol=1e-9,xtol=1e-9,
                gtol=1e-7,diff_step=1e-4)
            retry=fitter.solve(opt.x,pseudo)
            diagnostics.append(dict(cost=retry['cost'],nfev=opt.nfev,success=bool(opt.success),optimality=float(opt.optimality)))
            if retry['cost']<f['cost']:f=retry
            retried=True
        successful=any(d['success'] for d in diagnostics)
        row=dict(gene=gene,draw=b+1,SSE=f['cost'],exceeds_observed=f['cost']>=central['cost'],
            optimizer_success=successful,successful_starts=sum(d['success'] for d in diagnostics),
            retry= retried,total_nfev=sum(d['nfev'] for d in diagnostics))
        row.update({k:v for k,v in f['parameters'].items() if k!='pulses'})
        rows.append(row)
        allfits.append(dict(draw=b+1,parameters=f['parameters'],amplitudes=f['amplitudes'],
            offset=f['offset'],cost=f['cost'],optimization=diagnostics))
        if (b+1)%20==0 or b==B-1:
            print('null',gene,b+1,'/',B,'seconds',round(time.time()-t0,1),
                'optimizer_failures',sum(not r['optimizer_success'] for r in rows),flush=True)
            pd.DataFrame(rows).to_csv(OUT/f'null_{gene}.csv',index=False)
    count=sum(r['exceeds_observed'] for r in rows)
    result=dict(gene=gene,B=B,seed=2026091300+GENES.index(gene),dt=fitter.dt,
        observed_SSE=central['cost'],null_SSE_mean=np.mean([r['SSE'] for r in rows]),
        null_SSE_median=np.median([r['SSE'] for r in rows]),exceedances=count,
        conditional_lack_of_fit_p=(1+count)/(B+1),
        monte_carlo_se=np.sqrt(((1+count)/(B+1))*(1-(1+count)/(B+1))/(B+1)),
        optimizer_failures=sum(not r['optimizer_success'] for r in rows),
        retries=sum(r['retry'] for r in rows),starts=starts,scale=scale,
        elapsed_seconds=time.time()-t0,source_file=str(path.relative_to(ROOT)),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        interpretation='Conditional approximate descriptive lack-of-fit test; not cross-validation, proof of no overfitting, or biochemical parameter inference.',
        method='Within-stage paired residual bootstrap centered at the fitted bivariate mean; finite-n residual correction sqrt(n/(n-1)); fixed original channel scaling; both paired channels resampled together; same full mapped constrained cosine pulse model refitted.')
    (OUT/f'assessment_{gene}.json').write_text(json.dumps(serialize(result),indent=2))
    (OUT/f'null_fits_{gene}.json').write_text(json.dumps(serialize(allfits),indent=2))
    print('NULL COMPLETE',gene,'p',result['conditional_lack_of_fit_p'],flush=True)

def optimizer_check(gene):
    """Exactly two prespecified global-search checks; retain any better solutions."""
    fitter,fits,pairs,scale,y,path=source(gene)
    central=fitter.solve(fitter.pack(fits['mapped']['parameters']),y)
    result=json.loads((OUT/f'assessment_{gene}.json').read_text())
    allfits=json.loads((OUT/f'null_fits_{gene}.json').read_text())
    df=pd.read_csv(OUT/f'null_{gene}.csv')
    old_p=result['conditional_lack_of_fit_p'];checks=[];t0=time.time()
    rng=np.random.default_rng(result['seed'])
    residuals=[(p-p.mean(0))/scale*np.sqrt(len(p)/(len(p)-1)) for p in pairs]
    for b in range(2):
        pseudo=central['pred']+np.array([r[rng.integers(0,len(r),size=len(r))].mean(0) for r in residuals])
        before=float(df.loc[b,'SSE'])
        original=allfits[b]
        reconstructed=fitter.solve(fitter.pack(original['parameters']),pseudo)
        if abs(before-reconstructed['cost'])>1e-8:
            raise ValueError('Null draw reconstruction failed.')
        f=fitter.fit(pseudo,starts=[original['parameters']]+result['starts'],
            global_search=True,seed=2026091400+10*GENES.index(gene)+b)
        improved=f['cost']<before
        checks.append(dict(gene=gene,draw=b+1,local_SSE=before,
            global_check_SSE=f['cost'],improvement=max(0.,before-f['cost']),
            retained_improvement=improved,global_seed=2026091400+10*GENES.index(gene)+b))
        if improved:
            allfits[b]=dict(draw=b+1,parameters=f['parameters'],amplitudes=f['amplitudes'],
                offset=f['offset'],cost=f['cost'],optimization=f['optimization'],global_checked=True)
            df.loc[b,'SSE']=f['cost'];df.loc[b,'exceeds_observed']=f['cost']>=result['observed_SSE']
            for k in ['T','splice','decay','gain']:df.loc[b,k]=f['parameters'][k]
        print('global check',gene,b+1,'local',before,'global',f['cost'],
            'seconds',round(time.time()-t0,1),flush=True)
    count=int(df.exceeds_observed.sum());B=result['B']
    result.update(exceedances=count,conditional_lack_of_fit_p=(1+count)/(B+1),
        null_SSE_mean=float(df.SSE.mean()),null_SSE_median=float(df.SSE.median()),
        optimizer_spotchecks=2,p_before_spotchecks=old_p,
        optimizer_spotcheck_max_improvement=max(c['improvement'] for c in checks),
        optimizer_spotcheck_elapsed_seconds=time.time()-t0)
    result['monte_carlo_se']=np.sqrt(result['conditional_lack_of_fit_p']*(1-result['conditional_lack_of_fit_p'])/(B+1))
    result['optimizer_spotcheck_changed_p']=result['conditional_lack_of_fit_p']!=old_p
    df.to_csv(OUT/f'null_{gene}.csv',index=False)
    pd.DataFrame(checks).to_csv(OUT/f'optimizer_spotcheck_{gene}.csv',index=False)
    (OUT/f'assessment_{gene}.json').write_text(json.dumps(serialize(result),indent=2))
    (OUT/f'null_fits_{gene}.json').write_text(json.dumps(serialize(allfits),indent=2))
    print('CHECK COMPLETE',gene,'p',result['conditional_lack_of_fit_p'],flush=True)

def summary():
    rows=[]
    for gene in GENES:
        result=json.loads((OUT/f'assessment_{gene}.json').read_text())
        keys=['gene','B','dt','observed_SSE','null_SSE_median','null_SSE_mean',
            'exceedances','conditional_lack_of_fit_p','monte_carlo_se',
            'optimizer_failures','retries','optimizer_spotchecks',
            'optimizer_spotcheck_max_improvement','optimizer_spotcheck_changed_p']
        rows.append({k:result.get(k) for k in keys})
    pd.DataFrame(rows).to_csv(OUT/'assessment_summary.csv',index=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--gene',choices=GENES)
    ap.add_argument('--task',choices=['profile','null','check','summary'],required=True)
    ap.add_argument('--dt',type=float,default=None);ap.add_argument('--draws',type=int,default=199)
    args=ap.parse_args();OUT.mkdir(exist_ok=True)
    if args.task!='summary' and args.gene is None:ap.error('--gene is required for this task')
    if args.task=='summary':summary()
    elif args.task=='profile':profile(args.gene,args.dt)
    elif args.task=='null':null(args.gene,args.draws,args.dt)
    else:optimizer_check(args.gene)
