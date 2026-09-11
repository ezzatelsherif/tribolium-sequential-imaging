"""Independent adaptive-quadrature check of fitted pulse convolutions."""
import json
import numpy as np
import pandas as pd
from scipy.integrate import quad
from pulse_fit import ROOT,GENES,PulseFit
from pulse_forward import response,pulse_values
from position_kernels import cumulative,age_kernel

def reference(g,p,t,shape,mode):
    T=p['T'];sp=p['splice'];de=p['decay'];ts=T*g['splice_site_bp']/g['transcript_end_bp']
    end=t-(-12 if shape=='cosine' else -120)
    events=[0.,end,T,ts]
    events += [T*x/g['transcript_end_bp'] for ranges in [g['exon_targets'],g['intron_targets']] for interval in ranges for x in interval]
    events += [ts+k*sp for k in [.25,1.,4.,16.,64.]]+[T+k*de for k in [.25,1.,4.,16.,64.]]
    for q in p['pulses']:
        events += [t-q['center'],t-q['center']+q['rise'],t-q['center']-q['fall']]
    events=sorted(set(x for x in events if 0<=x<=end))
    ai,_,ae,_=cumulative(np.array(1e5),g,T,sp,de,mode)
    results=np.zeros((2,len(p['pulses'])))
    for j,q in enumerate(p['pulses']):
        for ch,area in enumerate([ai,ae]):
            fun=lambda age:float(age_kernel(age,g,T,sp,de,mode)[ch]*pulse_values(t-age,q['center'],q['rise'],q['fall'],shape))/area
            results[ch,j]=sum(quad(fun,a,b,epsabs=2e-11,epsrel=2e-10,limit=150)[0] for a,b in zip(events[:-1],events[1:]))
    return results

def main():
    rows=[]
    for gene in GENES:
        for shape in ['cosine','gaussian']:
            fits=json.loads((ROOT/f'fitresults/fit_{gene}_{shape}.json').read_text())
            for mode,f in fits.items():
                p=f['parameters'];g=PulseFit(gene).g
                for t in [4.,8.,11.,17.,20.,24.]:
                    exact=reference(g,p,t,shape,mode)
                    for dt in [.025,.0125,.00625]:
                        zi,ze=response(g,p,[t],shape,mode,dt);got=np.vstack([zi,ze])
                        rows.append(dict(gene=gene,shape=shape,mode=mode,stage_index=t,dt=dt,maximum_unit_amplitude_error=float(np.max(np.abs(got-exact)))))
    table=pd.DataFrame(rows);table.to_csv(ROOT/'fitresults/quadrature_validation.csv',index=False)
    print(table.groupby('dt').maximum_unit_amplitude_error.max().to_string())

if __name__=='__main__':main()
