"""Referee data checks (finite): closed forms vs recorded QSI data.
 (1) rho(delta)/rho(0), delta=1..4, all keys d=4..18 (binteg_d data)  [extends the author's d<=11, 72 cells]
 (2) Newton coefficients e(delta,s) from tower_k data (k<=14) vs closed forms F, all keys, incl. diagonal e-hat.
"""
import sys, os, time
ORIG = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, ORIG)
sys.path.insert(1, '.')
os.chdir('.')
import towlib as TL
from towlib import QSI, ZERO, keys_of_degree
import t1_direct as H
from fractions import Fraction as Fr
from math import comb

def toG(q):
    assert q.p1 == 0 and q.p3 == 0, q
    return H.G(q.p0, q.p2)

# (1) rho ratios
bad = 0; n = 0; t0 = time.time()
for dd in range(4, 19):
    for key in keys_of_degree(dd):
        p, q = key
        r0 = TL.rho(key, 0)
        for delta in range(1, 5):
            if delta > 2*dd-4: continue
            rr = TL.rho(key, delta) / r0
            n += 1
            if not (toG(rr) == H.arho(delta, dd, p)):
                bad += 1; print('rho mismatch', key, delta)
print(f'(1) rho(delta)/rho(0) closed form vs binteg data: d=4..18 all keys, delta=1..4: {n} cells, {bad} bad ({round(time.time()-t0,1)}s)')

# (1b) rho(0) = e00 closed form
bad = 0; n = 0
for dd in range(2, 19):
    for key in keys_of_degree(dd):
        p, q = key
        r0 = TL.rho(key, 0)
        # r0 = e00 * 1/(2 sqrt2)  -> r0 * 2sqrt2 should be e00 (Q(i))
        s = r0 * QSI(0, 2, 0, 0)    # 2*sqrt2
        n += 1
        if not (toG(s) == H.e00(p, q)):
            bad += 1; print('e00 mismatch', key)
print(f'(1b) rho(0) = K_d(-1)^p(1+[p>q])binom(d,p): d=2..18: {n} keys, {bad} bad')

# (2) Newton coefficients from tower data
def newton_from_data(key, delta, s):
    p, q = key; dd = p+q
    base = (2*p+1) if p == q else dd
    return sum((TL.u(p, q, base+t, delta).scale(Fr((-1)**(s-t)*comb(s, t))) for t in range(0, s+1)), ZERO), base

bad = 0; n = 0
for dd in range(4, 15):
    for key in keys_of_degree(dd):
        p, q = key
        for (delta, s) in [(2,1),(3,1),(4,1),(4,2)]:
            if delta > 2*dd-4: continue
            base = (2*p+1) if p == q else dd
            if base + s > 14: continue
            ed, base = newton_from_data(key, delta, s)
            if p > q:
                pred = H.e_abs(p, q, delta, s)
            else:
                pred = H.ehat(p, delta, s)
            n += 1
            # data values carry the factor 1/(2 sqrt2): multiply by 2sqrt2
            got = toG(ed * QSI(0, 2, 0, 0))
            if not (got == pred):
                bad += 1; print('newton mismatch', key, delta, s)
print(f'(2) Newton e(delta,s) (s>=1) from tower data vs closed forms (diag: e-hat): {n} cells, {bad} bad')
