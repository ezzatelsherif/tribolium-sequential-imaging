"""Fast pulse convolution with exact integrals of mapped biological age kernels.

Only the input pulse is piecewise-linearly approximated. Kernel transitions and
waiting times, including times smaller than dt, are integrated analytically.
"""
import numpy as np
from scipy.signal import fftconvolve
from position_kernels import cumulative

def pulse_values(t,c,rise,fall,shape='cosine'):
    t=np.asarray(t); z=t-c; w=np.where(z<0,rise,fall)
    if shape=='gaussian': return np.exp(-0.5*(z/w)**2)
    if shape!='cosine': raise ValueError(shape)
    return np.where(np.abs(z)<=w,np.cos(np.pi*z/(2*w))**2,0.)

def kernel_weights(g,T,splice,decay,mode,dt,span):
    a=dt*np.arange(int(np.ceil(span/dt))+2)
    f0i,f1i,f0e,f1e=cumulative(a,g,T,splice,decay,mode)
    weights=[]
    for f0,f1 in [(f0i,f1i),(f0e,f1e)]:
        d0=np.diff(f0); d1=np.diff(f1)
        left=(a[1:]*d0-d1)/dt
        right=(d1-a[:-1]*d0)/dt
        w=np.r_[left[0],right[:-1]+left[1:],right[-1]]
        # Roundoff in remote-tail differences can be a few 1e-10.
        weights.append(np.maximum(w,0.))
    ai,_,ae,_=cumulative(np.array(1e5),g,T,splice,decay,mode)
    return weights[0]/ai,weights[1]/ae

def response(g,p,at,shape='cosine',mode='mapped',dt=.05):
    at=np.asarray(at,float)
    start=-12. if shape=='cosine' else -120.
    end=max(26.,float(np.max(at)))
    grid=start+dt*np.arange(int(np.ceil((end-start)/dt))+1)
    wi,we=kernel_weights(g,p.get('T',0.),p['splice'],p['decay'],mode,dt,end-start)
    outi=[];oute=[]
    for q in p['pulses']:
        u=pulse_values(grid,q['center'],q['rise'],q['fall'],shape)
        yi=fftconvolve(u,wi)[:len(grid)];ye=fftconvolve(u,we)[:len(grid)]
        outi.append(np.interp(at,grid,yi));oute.append(np.interp(at,grid,ye))
    return np.column_stack(outi),np.column_stack(oute)
