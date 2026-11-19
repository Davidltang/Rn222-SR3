"""Toy feasibility test: can Rn-222 be counted with AFT-only correction (no r2 correction)?
Takes fully corrected S1-only areas (perbin_export/areas_by_partition.npz), UNDOES Correction 3 with an
approximate radial light-collection model g(r2) (~9% loss at the wall, from the Corr-3 fit output:
median area 71.9k at r2=0 -> 65.6k), then compares rate-extraction methods against the full-correction answer."""
import numpy as np, sys
from scipy.optimize import curve_fit
z=np.load('perbin_export/areas_by_partition.npz')
G=lambda x,A,m,s: A*np.exp(-0.5*((x-m)/s)**2)
def model(x,Aw,mw,sw,A1,m1,s1,A2,m2,s2,b):  # wall + Rn222 + Po218 + flat
    return G(x,Aw,mw,sw)+G(x,A1,m1,s1)+G(x,A2,m2,s2)+b
def fit_full(a,bw=250,lo=38e3,hi=62e3,mu=50.5e3):
    e=np.arange(lo,hi+bw,bw); h,_=np.histogram(a,e); x=0.5*(e[1:]+e[:-1])
    p0=[h.max()*.1,0.92*mu,1.5e3,h.max(),mu,0.9e3,h.max()*.9,mu*6.0/5.49,0.9e3,1]
    lb=[0,0.85*mu,300,0,0.97*mu,300,0,1.05*mu,300,0]; ub=[np.inf,0.96*mu,4e3,np.inf,1.03*mu,4e3,np.inf,1.15*mu,4e3,np.inf]
    p,_=curve_fit(model,x,h,p0=p0,bounds=(lb,ub),sigma=np.sqrt(h+1),maxfev=40000)
    return p,bw
def n_of(A,s,bw): return A*s*np.sqrt(2*np.pi)/bw
def naive_window(a,p):  # Caio-style: count +-2 sigma around Rn222 fit, no leakage model for wall
    m1,s1=p[4],p[5]; return ((a>m1-2*s1)&(a<m1+2*s1)).sum()/0.9545
def upper_side(a,p,k=1.5):
    m1,s1=p[4],p[5]; from scipy.stats import norm; f=norm.cdf(k)-0.5
    return ((a>m1)&(a<m1+k*s1)).sum()/f
def po218_route(a,p,k=1.5,ratio=0.93):
    m2,s2=p[7],p[8]; from scipy.stats import norm; f=norm.cdf(k)-0.5
    return ((a>m2)&(a<m2+k*s2)).sum()/f/ratio
for wall_loss,pw in [(0.088,1.0),(0.088,2.0)]:
    print(f"\n=== toy g(r2)=1-{wall_loss}*(r2/50)^{pw} ===")
    print(f"{'part':4s} {'truth(full)':>11s} {'naive+-2s':>10s} {'4G fit':>8s} {'upper1.5s':>10s} {'Po218/0.93':>10s}  sigma_full->aftonly")
    for P in ['P1','P2','P3','P4']:
        a=z[P]; r2=np.clip(z[P+'_r2'],0,50)
        pt,bw=fit_full(a); T=n_of(pt[3],pt[5],bw)
        aa=a*(1-wall_loss*(r2/50)**pw)
        mu=pt[4]*np.median(1-wall_loss*(r2/50)**pw)
        try: pa,bw=fit_full(aa,mu=mu)
        except Exception as ex: print(P,'fit fail',ex); continue
        res=[naive_window(aa,pa),n_of(pa[3],pa[5],bw),upper_side(aa,pa),po218_route(aa,pa)]
        print(f"{P:4s} {T:11.0f} "+" ".join(f"{v:7.0f}({v/T-1:+.0%})" for v in res)+f"  {pt[5]:.0f}->{pa[5]:.0f}")
