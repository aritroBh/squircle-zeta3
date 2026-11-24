# Planner check of TOW1R: (a) exact kernels R20, R11 at new odd k with the paper's reducer;
# (b) the slope constant identified by PSLQ (not assumed) from the slope kernels alone.
import sys, time
sys.dont_write_bytecode = True
WTS = "<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0, WTS)
from fractions import Fraction as Fr
import mpmath as mp
from fastred import *
from buildintegrand import build_integrand_pf

def claimed(k):
    r20 = {ZERO: QSI(Fr(1,4), Fr(-1,4)), AQ: QSI(Fr(1,4)), -AQ: QSI(Fr(2*k-7,4)), -BQ: QSI(Fr(-(k-2),2))}
    r11 = {ZERO: QSI(Fr(2*k-5,4), Fr(1,4)), -AQ: QSI(Fr(-(2*k-5),2))}
    return r20, r11

def same(P, E):
    return all(P.get(r, ZERO) == c for r, c in E.items()) and all(c.is_zero() or r in E for r, c in P.items())

for k in [int(x) for x in sys.argv[1:]] or [23]:
    t0 = time.time()
    fin = reduce_k_fast(k, build_integrand_pf(k))
    e20, e11 = claimed(k)
    print(f"k={k}: R20 {same(fin[(2,0)][0], e20)}  R11 {same(fin[(1,1)][0], e11)}  higher {fin[(2,0)][1:]} {fin[(1,1)][1:]}  ({time.time()-t0:.0f}s)", flush=True)

mp.mp.dps = 110
s2 = mp.sqrt(2); a = 1 + s2; b = s2 - 1
t = lambda v: (v*v + 2*v - 1)/(4*v)
om = lambda v: (1 + 2*v - v*v)/(4*v)
S20 = lambda v: 1/(2*(v + a)) - 1/(2*(v + b))
S11 = lambda v: 1/(2*v) - 1/(v + a)
slope = mp.quad(lambda v: S20(v)*mp.log(t(v))**2 + S11(v)*mp.log(t(v))*mp.log(om(v)), [b, 1, a])
rel = mp.pslq([slope, mp.zeta(3), mp.pi**2*mp.log(2), mp.log(2)**3], maxcoeff=10**6, maxsteps=10**6)
print("slope =", mp.nstr(slope, 30))
print("PSLQ [slope, z3, pi^2 log2, log^3 2]:", rel, "(claim: slope = -(z3/8 + pi^2 log2/24 + log^3 2/12) -> [24, 3, 1, 2])")
