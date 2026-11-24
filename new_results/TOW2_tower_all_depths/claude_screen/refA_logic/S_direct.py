"""Referee: direct numeric check of S_m^{(l)}(t*)/S_m(t*) = N_l/Pi_l from the DEFINITION of S_m (genintegrand.S_poly),
for all d=3..DMAX and ALL m=0..d (incl. m<l boundary and diagonal), l=1..4. Exact Gaussian-rational arithmetic."""
import sys, time
ORIG = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, ORIG); sys.path.insert(1, '.')
from fractions import Fraction as Fr
import genintegrand as GI
import tow2lib as T
import t1_direct as H
G = H.G
tstar = G(Fr(1,2), Fr(1,2))
def deriv_at(P, l):
    # l-th derivative of poly P (list low->high Fr) at tstar, Gaussian rational
    # synthetic: Taylor coefficients via repeated Horner
    n = len(P)
    coeffs = [G(c) for c in P]
    res = []
    cur = coeffs
    fact = 1
    out = []
    for lev in range(l+1):
        # evaluate cur at tstar and replace cur by its derivative polynomial coefficients
        val = G(0)
        for c in reversed(cur): val = val*tstar + c
        out.append(val)
        cur = [cur[j]*j for j in range(1, len(cur))]
        if not cur: cur = [G(0)]
    return out
import sympy as sp
def Nval(l, dd, mm):
    e = T.N[l].subs({T.d: dd, T.m: mm})
    pi = T.Pi(l).subs({T.d: dd})
    return H.sp2G(e/pi)
DMAX = int(sys.argv[1]) if len(sys.argv)>1 else 14
bad = {1:0,2:0,3:0,4:0}; n=0; t0=time.time()
for dd in range(int(sys.argv[2]) if len(sys.argv)>2 else 5, DMAX+1):
    for mm in range(0, dd+1):
        P = GI.S_poly(dd, mm)
        vals = deriv_at(P, 4)
        S0 = vals[0]
        for l in (1,2,3,4):
            ratio = vals[l] / S0
            n += 1
            if not (ratio == Nval(l, dd, mm)):
                bad[l] += 1
                if bad[l] <= 3: print('mismatch l=',l,'d=',dd,'m=',mm)
print(f'S^(l)(t*)/S(t*) vs N_l/Pi_l from definition: d=5..{DMAX}, all m: {n} (d,m,l) cells; bad by l = {bad}  ({round(time.time()-t0,1)}s)')
