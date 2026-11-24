"""Referee re-implementation: ALL non-adjacent Newton identities N(delta,s), delta=0..4, as exact rational identities in (d,m).
Uses tow2lib.arho_simple (pickle) for rho and the candidate e-family; tau table from tow2lib (checked separately).
Ratios rA0, rB0 are re-derived from e00 = K_d (-1)^p (1+[p>q]) binom(d,p), K_{d+1}/K_d = -i(2d-1)/(d+1)."""
import sys, time
sys.path.insert(0,'.')
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
i = sp.I
# re-derive ratios: e00(p+1,q)/e00(p,q) with p=m, q=d-m:  binom(d+1,m+1)/binom(d,m) = (d+1)/(m+1); sign flips (-1); K ratio
Kr = -i*(2*d-1)/(d+1)
rA0 = sp.simplify(Kr * (-1) * (d+1)/(m+1))
rB0 = sp.simplify(Kr * (+1) * (d+1)/(d-m+1))      # binom(d+1,m)/binom(d,m) = (d+1)/(d-m+1), same sign
assert sp.simplify(rA0 - i*(2*d-1)/(m+1)) == 0 and sp.simplify(rB0 + i*(2*d-1)/(d-m+1)) == 0
print('rA0, rB0 re-derived OK')
def F(delta, s, dd, mm):
    if s > delta//2: return sp.Integer(0)
    if s == 0: return T.arho_simple(delta).subs({d: dd, m: mm}, simultaneous=True) if delta else sp.Integer(1)
    if (delta,s) == (2,1): return -1/(2*dd-3)
    if (delta,s) == (3,1): return i*(dd-3)/(2*dd-3)
    if (delta,s) == (4,1): return (dd**3 - 9*dd**2 + sp.Rational(113,4)*dd - i*dd/2 + i*mm - 29)/((2*dd-3)*(2*dd-5))
    if (delta,s) == (4,2): return sp.Rational(2)/((2*dd-3)*(2*dd-5))
    raise ValueError
def Lam(delta, s):
    tot = 0
    for dp in range(delta+1):
        tot += T.tau1f(delta-dp)*(m+1)/(2*d-dp-1) * F(dp, s, d+1, m+1) * rA0
        tot += T.tau2f(delta-dp)*(d-m+1)/(2*d-dp-1) * F(dp, s, d+1, m) * rB0
    return tot
t0=time.time()
for (delta, s) in [(0,0),(1,0),(2,0),(2,1),(3,0),(3,1),(4,0),(4,1),(4,2)]:
    lhs = F(delta, s, d, m) + F(delta, s+1, d, m)
    rem = sp.cancel(sp.together(lhs - Lam(delta, s)))
    print(f'NA N({delta},{s}) remainder = {rem}   ({round(time.time()-t0,1)}s)', flush=True)
