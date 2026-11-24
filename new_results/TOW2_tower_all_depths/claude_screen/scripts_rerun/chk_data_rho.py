# Independent check of rho(delta) closed forms vs recorded initial data (binteg pickles)
import sympy as sp, time
from towlib import rho, keys_of_degree
import tow2lib as T
d, m = T.d, T.m
t0 = time.time(); n = 0; bad = []
for dd in range(2, 16):
    for key in keys_of_degree(dd):
        e00 = rho(key, 0)
        for delta in range(1, min(4, 2*dd-4)+1):
            r = rho(key, delta)/e00
            val = T.arho_simple(delta).subs({d: dd, m: key[0]})
            val = sp.expand(val)
            re, im = sp.re(val), sp.im(val)
            ok = (r.p1 == 0 and r.p3 == 0 and sp.Rational(r.p0.numerator, r.p0.denominator) == re and sp.Rational(r.p2.numerator, r.p2.denominator) == im)
            n += 1
            if not ok: bad.append((key, delta, r.p0, r.p1, r.p2, r.p3, val))
print(f'rho(delta)/rho(0) data vs arho_simple: {n} cells, d=2..15, delta=1..min(4,2d-4), {len(bad)} bad ({time.time()-t0:.1f}s)')
print(bad[:5])
