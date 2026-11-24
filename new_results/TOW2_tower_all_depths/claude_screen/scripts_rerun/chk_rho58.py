# Bonus: does the Bell machinery of tow2lib (using proved N1..N4) reproduce recorded rho(delta) for delta=5..8 ?
import sympy as sp, time
from towlib import rho, keys_of_degree
import tow2lib as T
d, m = T.d, T.m; y = sp.Symbol('y')
Lsym = {r: T.L(r) for r in range(1, 9)}
t0 = time.time(); n = 0; bad = []
for dd in range(6, 13):
    for key in keys_of_degree(dd):
        Lnum = {r: sp.radsimp(sp.expand(Lsym[r].subs({d: dd, m: key[0]}))) for r in range(1, 9)}
        ser = sp.series(sp.exp(sum(Lnum[r]*y**r for r in range(1, 9))), y, 0, 9).removeO()
        e00 = rho(key, 0)
        for delta in range(5, min(8, 2*dd-4)+1):
            val = sp.expand(sp.radsimp(sp.expand(ser.coeff(y, delta))))
            r = rho(key, delta)/e00
            ok = (r.p1 == 0 and r.p3 == 0 and sp.Rational(r.p0.numerator, r.p0.denominator) == sp.re(val) and sp.Rational(r.p2.numerator, r.p2.denominator) == sp.im(val))
            n += 1
            if not ok: bad.append((key, delta))
print(f'rho(5..8) via Bell(L_1..L_8, proved N1..N4) vs recorded data: {n} cells, {len(bad)} bad  ({time.time()-t0:.1f}s)')
print(bad[:8])
