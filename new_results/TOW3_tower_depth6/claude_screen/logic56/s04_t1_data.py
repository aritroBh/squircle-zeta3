# s04: validate (i) e00 formula K_d(-1)^p(1+[p>q])binom(d,p), (ii) my T1 implementation (incl. fold, true binomials)
# against the exact TOWN1 data (independent of every TOW3 closed form).
import sys, time
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import sympy as sp
import tow2lib as T2
from t1lib import *
t0 = time.time()
def gq(z):
    z = sp.nsimplify(z); return GQ(Fr(int(sp.re(z).p), int(sp.re(z).q)), Fr(int(sp.im(z).p), int(sp.im(z).q)))
taus = Taus([gq(T2.tau1f(n)) for n in range(8)], [gq(T2.tau2f(n)) for n in range(8)])
H = load_data()
keys = [k for k in H if k != '_meta']
print('data keys', len(keys), 'max degree', max(p+q for p, q in keys), flush=True)
# (i) e00
bad = 0; n = 0
for (p, q) in keys:
    dd = p+q
    if dd < 2: continue
    if p == q: kk = 2*p
    else: kk = dd
    if not data_has(p, q, kk, 0): continue
    n += 1
    if not (data_get(p, q, kk, 0) == E00(p, q)): bad += 1
print(f'(i) rho(0) data == (sqrt2/4)*K_d(-1)^p(1+[p>q])C(d,p)*...: {n} keys, {bad} bad', flush=True)
# (ii) T1 vs data
tested = 0; bad = []; bydelta = {}
for (p, q) in sorted(keys, key=lambda x: (x[0]+x[1], x)):
    dd = p+q
    if dd < 3 or dd > 30: continue
    for k in range(dd+1, min(dd+7, 52)+1):
        for delta in range(0, min(6, 2*dd-4)+1):
            # need all parents in data
            need = [(p+1, q)] + ([(p, q+1)] if p > q else []) + ([(p, p)] if q == p-1 else [])
            if not all(data_has(a, b, k, 0) for a, b in need): continue
            if not data_has(p, q, k, delta): continue
            lhs = data_get(p, q, k, delta)
            rhs = T1rhs(data_get, taus, p, q, k, delta)
            tested += 1; bydelta[delta] = bydelta.get(delta, 0) + 1
            if not (lhs == rhs): bad.append((p, q, k, delta))
print(f'(ii) my T1 (fold w/ true binomials) vs TOWN1 data, d=3..30, k=d+1..d+7, delta<=6: {tested} cells, {len(bad)} bad', bydelta, flush=True)
if bad[:10]: print(bad[:10])
print('done', round(time.time()-t0, 1), 's')
