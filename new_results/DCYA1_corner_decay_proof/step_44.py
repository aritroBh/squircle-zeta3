import os
cwd='<project>/research/08_theory/DCYA1_corner_decay_proof/'
os.chdir(cwd)
module = '''"""Reliable Taylor-coefficient evaluator for the corner function M(u;alpha).
M(u;a) = int_0^1 [ (1+g)^(-1-u) sqrt(1+g^(2-2u)-2cos(a)g^(1-u)) - 1 ] g^(u-1) dg = sum_k M_k u^k.
Method: exact formal u-series of the integrand at each tanh-sinh node (shared across all k),
then tanh-sinh quadrature per coefficient. Handles the g^(u-1) endpoint singularity.
Set mp.mp.dps high enough (grows with k; k=100->~140, k=300->~350)."""
import mpmath as mp
import corner_series as cs

def mk_all(a, N, pts=None, h=None, n=None):
    """Return ([M_0..M_N], nnodes). Call with mp.mp.dps set."""
    if pts is None:
        pts=[mp.mpf(0),mp.mpf(1)/64,mp.mpf(1)/8,mp.mpf(1)/2,mp.mpf(1)]
    if h is None: h=mp.mpf('0.015')
    if n is None: n=350
    c=mp.cos(a); cache={}
    def Fk(g,k):
        if g not in cache: cache[g]=cs.integrand_series(g,c,N)
        return cache[g][k]
    tiny=mp.mpf(10)**(-mp.mp.dps+10)
    def ts_integrate(f,a_,b_):
        half=(b_-a_)/2; mid=(a_+b_)/2; S=mp.mpf(0)
        for i in range(-n,n+1):
            t=i*h; s=mp.sinh(t); th=mp.tanh(mp.pi/2*s)
            x=mid+half*th
            if x<=a_+tiny or x>=b_-tiny: continue
            dx=half*(mp.pi/2)*mp.cosh(t)/(mp.cosh(mp.pi/2*s)**2)
            S+=f(x)*dx
        return S*h
    M=[mp.mpf(0)]*(N+1)
    for k in range(N+1):
        M[k]=mp.fsum(ts_integrate(lambda g: Fk(g,k), pts[i],pts[i+1]) for i in range(len(pts)-1))
    return M,len(cache)

def lead(a,k):
    N=k+1; A=mp.sqrt(2*a*mp.sin(a)); th=mp.sqrt(2*a*N)
    return -(A/N)*mp.e**(-th)*mp.cos(th), A, th

def resid_ratio(Mk,k,a):
    lv,A,th=lead(a,k); N=k+1
    sc=N**(-mp.mpf(3)/2)*mp.e**(-th)
    return abs(Mk-lv)/sc
'''
open(cwd+'mk_eval.py','w').write(module)
# log progress
with open(cwd+'progress.log','a') as f:
    f.write("step: built reliable M_k evaluator (custom series+tanh-sinh); confirmed true M_100(pi/2)=-1.7723168073590541001e-10 at dps=140 (matches K2). FFT method has noise floor, abandoned.\n")
print("mk_eval.py written; progress logged.")
print(open(cwd+'progress.log').read())
