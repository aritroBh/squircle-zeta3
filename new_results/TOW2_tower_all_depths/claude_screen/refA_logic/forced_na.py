"""Fit-free derivation (NA keys): e(delta,s+1) := Lambda^(delta)(s) - e(delta,s), e(delta,0)=rho(delta), recursion in the same (d,m)-family.
Shows the 'fitted' closed forms are FORCED, and truncation e(delta, floor(delta/2)+1) = 0 is exactly the top identity."""
import sys, time
sys.path.insert(0,'.')
import sympy as sp
import tow2lib as T
d, m, i = T.d, T.m, sp.I
Kr = -i*(2*d-1)/(d+1)
rA0 = i*(2*d-1)/(m+1); rB0 = -i*(2*d-1)/(d-m+1)
from functools import lru_cache
@lru_cache(None)
def E(delta, s, dd=None, mm=None):
    """forced e(delta,s)/e00 at (dd,mm) (default symbolic d,m), as rational function; computed WITHOUT the fitted closed forms."""
    dd = d if dd is None else dd; mm = m if mm is None else mm
    if s > delta//2 + 1: return sp.Integer(0)  # only used to expose the truncation check; real zero beyond is imposed in caller
    if s == 0:
        return sp.Integer(1) if delta == 0 else T.arho_simple(delta).subs({d: dd, m: mm}, simultaneous=True)
    # forced: e(delta,s) = Lambda(delta,s-1) - e(delta,s-1)
    return sp.cancel(sp.together(Lam(delta, s-1, dd, mm) - E(delta, s-1, dd, mm)))
def Lam(delta, s, dd, mm):
    tot = 0
    rA = sp.Symbol('rA'); 
    rA0_ = (i*(2*dd-1)/(mm+1)); rB0_ = (-i*(2*dd-1)/(dd-mm+1))
    for dp in range(delta+1):
        if s > dp//2: continue   # parent's e(dp,s)=0 for s>floor(dp/2) (truncation inductively established for dp<delta; for dp=delta top is the check)
        tot += T.tau1f(delta-dp)*(mm+1)/(2*dd-dp-1) * E(dp, s, dd+1, mm+1) * rA0_
        tot += T.tau2f(delta-dp)*(dd-mm+1)/(2*dd-dp-1) * E(dp, s, dd+1, mm) * rB0_
    return tot
claimed = {(2,1): -1/(2*d-3), (3,1): i*(d-3)/(2*d-3),
           (4,1): (d**3 - 9*d**2 + sp.Rational(113,4)*d - i*d/2 + i*m - 29)/((2*d-3)*(2*d-5)),
           (4,2): sp.Rational(2)/((2*d-3)*(2*d-5))}
t0=time.time()
for key,val in claimed.items():
    forced = E(*key)
    print(f'forced e{key} - claimed = ', sp.cancel(sp.together(forced - val)), f' ({round(time.time()-t0,1)}s)', flush=True)
# truncation checks: next-level forced coefficient must vanish
for (delta, s) in [(2,2),(3,2),(4,3)]:
    print(f'forced e({delta},{s}) (must be 0) = ', sp.cancel(sp.together(E(delta, s))), f' ({round(time.time()-t0,1)}s)', flush=True)
