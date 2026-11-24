"""step_14 (TOW5R): rho(9), rho(10) = a_9, a_10 vs data — own Q(sqrt2,i) division.
 (a) TOW5's data_rho_D10_d7_40.pkl: ALL cells (not a sample), ratio asserted in Q(i).
 (b) WTS ground truth binteg_d.pkl (d = 7..18) via towlib.rho — the actual reduction data,
     independent of TOW5's engine. rho(delta)/rho(0) == a_delta.
"""
import os, sys, pickle
from fractions import Fraction as Fr
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE + '/tow5_copy')
_WTS = '<project>/research/08_theory/WTS_weight_k_structure/theory'
sys.path.insert(0, _WTS)
import sympy as sp
import tow5lib as T5
d, m = T5.d, T5.m

# ---- own Q(sqrt2,i) division ----
def gm(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gd(a, b):
    den = b[0]**2+b[1]**2; return ((a[0]*b[0]+a[1]*b[1])/den, (a[1]*b[0]-a[0]*b[1])/den)
def qsi_div(x, y):
    A = (x[0], x[2]); B = (x[1], x[3]); C = (y[0], y[2]); D = (y[1], y[3])
    C2 = (gm(C, C)[0]-2*gm(D, D)[0], gm(C, C)[1]-2*gm(D, D)[1])
    nA = (gm(A, C)[0]-2*gm(B, D)[0], gm(A, C)[1]-2*gm(B, D)[1])
    nB = (gm(B, C)[0]-gm(A, D)[0], gm(B, C)[1]-gm(A, D)[1])
    return (gd(nA, C2), gd(nB, C2))

# (a) TOW5 data, all cells
RHO = pickle.load(open(_HERE + '/tow5_copy/data_rho_D10_d7_40.pkl', 'rb'))
res = {}
for delta in [9, 10]:
    aa = T5.arho(delta)
    bad = 0; n = 0
    for dd in RHO:
        if 2*dd-4 < delta: continue
        for key in RHO[dd]:
            p_, q_ = key
            rv = qsi_div(RHO[dd][key][delta], RHO[dd][key][0])
            cv = sp.expand(aa.subs({d: dd, m: p_}))
            exp = (Fr(sp.re(cv)), Fr(sp.im(cv)))
            n += 1
            if rv[1] != (Fr(0), Fr(0)) or rv[0] != exp: bad += 1
    res[delta] = (n, bad)
    print(f'(a) TOW5 data rho({delta})/rho(0) == a_{delta}: {n} cells, bad = {bad}', flush=True)

# (b) WTS binteg ground truth via towlib.rho
sys.path.insert(0, '<project>/research/08_theory/TOW2_tower_all_depths')
import towlib
resb = {}
for delta in [9, 10]:
    aa = T5.arho(delta)
    bad = 0; n = 0
    for dd in range(7, 19):
        if 2*dd-4 < delta: continue
        for key in towlib.keys_of_degree(dd):
            r0 = towlib.rho(key, 0); rd = towlib.rho(key, delta)
            rv = qsi_div((rd.p0, rd.p1, rd.p2, rd.p3), (r0.p0, r0.p1, r0.p2, r0.p3))
            cv = sp.expand(aa.subs({d: dd, m: key[0]}))
            exp = (Fr(sp.re(cv)), Fr(sp.im(cv)))
            n += 1
            if rv[1] != (Fr(0), Fr(0)) or rv[0] != exp: bad += 1
    resb[delta] = (n, bad)
    print(f'(b) WTS binteg rho({delta})/rho(0) == a_{delta}: {n} cells, bad = {bad}', flush=True)
ok = all(b == 0 for _, b in res.values()) and all(b == 0 for _, b in resb.values())
print('STEP_14 VERDICT:', 'PASS' if ok else 'FAIL', flush=True)
assert ok
