"""Constrained one/two-pulse joint fits. See fitdesign.md for fixed choices."""
from pathlib import Path
import argparse,json,time
import numpy as np
import pandas as pd
from scipy.optimize import nnls,least_squares,differential_evolution
from pulse_forward import response,pulse_values

ROOT=Path(__file__).resolve().parent
X=np.arange(1.,25.); GENES=['Kr','mlpt','svb']
BOUNDS={'T':(.005,6.),'splice':(.005,6.),'decay':(.01,6.),'gain':(.2,5.)}

def load():
    df=pd.read_csv(ROOT/'data/background_corrected_pairs.csv');df=df[df.complete_pair]
    return {g:[p[['intron','exon']].to_numpy() for _,p in df[df.gene==g.upper()].groupby('stage_index')] for g in GENES}

def geometry(gene):return json.loads((ROOT/'data/geometry.json').read_text())['genes'][gene]

def serialize(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,(np.floating,np.integer)):return value.item()
    if isinstance(value,dict):return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(v) for v in value]
    return value

class PulseFit:
    def __init__(self,gene,mode='mapped',shape='cosine',dt=.025,geom=None):
        self.gene=gene;self.mode=mode;self.shape=shape;self.dt=dt
        self.g=geometry(gene) if geom is None else geom
        self.npulse=1 if gene=='Kr' else 2
        self.names=(['T'] if mode=='mapped' else [])+['splice','decay','gain']
        self.bounds=[tuple(np.log(BOUNDS[k])) for k in self.names]
        for j in range(self.npulse):
            self.names += [f'center{j+1}',f'rise{j+1}',f'fall{j+1}']
            self.bounds += [(0.,12.) if j==0 else (12.,26.),(.25,12.),(.25,12.)]
    def unpack(self,theta):
        flat={k:float(np.exp(v) if k in BOUNDS else v) for k,v in zip(self.names,theta)}
        p={k:flat[k] for k in BOUNDS if k in flat};p.setdefault('T',0.)
        p['pulses']=[{k:flat[f'{k}{j+1}'] for k in ['center','rise','fall']} for j in range(self.npulse)]
        return p
    def pack(self,p):
        flat={k:v for k,v in p.items() if k!='pulses'}
        flat.setdefault('T',.005)
        for j,q in enumerate(p['pulses']):flat.update({f'{k}{j+1}':v for k,v in q.items()})
        x=np.array([np.log(max(flat[k],BOUNDS[k][0])) if k in BOUNDS else flat[k] for k in self.names])
        return np.clip(x,np.array(self.bounds)[:,0],np.array(self.bounds)[:,1])
    def matrices(self,p,at=X):
        ki,ke=response(self.g,p,at,self.shape,self.mode,self.dt)
        A=np.stack([ki,ke*p['gain']],axis=1).reshape(-1,self.npulse)
        return np.column_stack([A,np.tile(np.eye(2),(len(at),1))])
    def solve(self,theta,y,at=X):
        p=self.unpack(theta);A=self.matrices(p,at)
        norms=np.linalg.norm(A,axis=0);valid=norms>1e-12
        coef=np.zeros(A.shape[1]);coef[valid]=nnls(A[:,valid]/norms[valid],np.asarray(y).ravel(),maxiter=500)[0]/norms[valid]
        pred=(A@coef).reshape(-1,2);residual=pred-y
        return dict(cost=float(np.sum(residual**2)),parameters=p,theta=np.asarray(theta),amplitudes=coef[:self.npulse],offset=coef[-2:],pred=pred,residual=residual,mode=self.mode,shape=self.shape,gene=self.gene,dt=self.dt)
    def predict(self,fit,at=X):
        coef=np.r_[fit['amplitudes'],fit['offset']]
        return (self.matrices(fit['parameters'],at)@coef).reshape(-1,2)
    def fit(self,y,starts=None,global_search=True,seed=10,maxiter=None,at=X):
        y=np.asarray(y);lo,hi=np.array(self.bounds).T
        objective=lambda x:self.solve(x,y,at)['cost']
        residual=lambda x:self.solve(x,y,at)['residual'].ravel()
        candidates=[]
        if global_search:
            de=differential_evolution(objective,self.bounds,seed=seed,popsize=8,maxiter=200 if maxiter is None else maxiter,tol=1e-7,polish=False,updating='immediate')
            candidates.append((de.fun,de.x))
            for j in np.argsort(de.population_energies)[:3]:candidates.append((de.population_energies[j],de.population[j]))
        for s in starts or []:
            x=self.pack(s) if isinstance(s,dict) else np.clip(s,lo,hi)
            candidates.append((objective(x),x))
        if not candidates:
            default=dict(T=.3,splice=.2,decay=.5,gain=1.1,pulses=[dict(center=5.,rise=4.,fall=5.)])
            if self.npulse==2:default['pulses'].append(dict(center=17.,rise=5.,fall=5.))
            x=self.pack(default);candidates=[(objective(x),x)]
        initial=list(sorted(candidates,key=lambda q:q[0]))
        # Optimize all starts, because biochemical/input trade-offs cause multiple basins.
        diagnostics=[]
        for _,x in initial:
            opt=least_squares(residual,np.clip(x,lo,hi),bounds=(lo,hi),max_nfev=300 if global_search else 100,ftol=1e-10 if global_search else 1e-7,xtol=1e-9,gtol=1e-8 if global_search else 1e-6,diff_step=1e-4)
            cost=float(np.sum(opt.fun**2));candidates.append((cost,opt.x))
            diagnostics.append(dict(cost=cost,nfev=opt.nfev,success=bool(opt.success),optimality=float(opt.optimality)))
        best=min(candidates,key=lambda q:q[0]);out=self.solve(best[1],y,at)
        out['optimization']=diagnostics
        out['bound_hits']=[k for k,v,l,h in zip(self.names,out['theta'],lo,hi) if min(v-l,h-v)<.001*(h-l)]
        return out

def main(gene,shape='cosine',use_saved_starts=True):
    out=ROOT/'fitresults';out.mkdir(exist_ok=True)
    pairs=load()[gene];means=np.array([p.mean(0) for p in pairs]);scale=means.max(0);y=means/scale
    results={};t0=time.time()
    for mode in ['collapsed','mapped']:
        fitter=PulseFit(gene,mode,shape,dt=.05)
        starts=[]
        prior=out/f'fit_{gene}_{shape}.json'
        if use_saved_starts and prior.exists():
            saved=json.loads(prior.read_text())
            if mode in saved:starts.append(saved[mode]['parameters'])
        if mode=='mapped':
            for t in [.005,.05,.3,1.,3.]:starts.append(dict(results['collapsed']['parameters'],T=t))
        coarse=fitter.fit(y,starts=starts,seed=20260911+GENES.index(gene)+100*(shape=='gaussian'))
        fitter=PulseFit(gene,mode,shape,dt=.025)
        f=fitter.fit(y,starts=[coarse['parameters']],global_search=False)
        f['coarse_mesh_initial_cost']=coarse['cost']
        f['coarse_optimization']=coarse['optimization']
        f['scale']=scale;f['y']=y;f['elapsed_seconds']=time.time()-t0
        results[mode]=f
        (out/f'fit_{gene}_{shape}.json').write_text(json.dumps(serialize(results),indent=2))
        print(gene,shape,mode,'seconds',round(time.time()-t0),'SSE',f['cost'],'params',f['parameters'],'bounds',f['bound_hits'],flush=True)
    rows=[]
    for mode,f in results.items():
        for i in range(24):rows.append(dict(gene=gene,shape=shape,mode=mode,stage_index=i+1,intron=means[i,0],exon=means[i,1],intron_fit=f['pred'][i,0]*scale[0],exon_fit=f['pred'][i,1]*scale[1]))
    pd.DataFrame(rows).to_csv(out/f'curves_{gene}_{shape}.csv',index=False)
    return results

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--gene',choices=GENES,required=True);ap.add_argument('--shape',choices=['cosine','gaussian'],default='cosine')
    ap.add_argument('--fresh',action='store_true',help='Ignore saved fit parameters as warm starts; retain fixed random seeds.')
    args=ap.parse_args();main(args.gene,args.shape,use_saved_starts=not args.fresh)
