# Bonus (fast numeric version): rho(5..8)/rho(0) from the Bell machinery with the proved N1..N4 vs recorded data
import sympy as sp, time
from towlib import rho, keys_of_degree
import tow2lib as T
i = sp.I; y = sp.Symbol('y')
def Lnum(r, dd, mm):
    R = {l: (T.N[l]/T.Pi(l)).subs({T.d: dd, T.m: mm}) for l in range(1, 5)}
    R = {l: sp.nsimplify(v) if False else sp.simplify(v) for l, v in R.items()}
    sl = {1: R[1], 2: R[2]-R[1]**2, 3: R[3]-3*R[1]*R[2]+2*R[1]**3, 4: R[4]-4*R[1]*R[3]-3*R[2]**2+12*R[1]**2*R[2]-6*R[1]**4}
    p12 = (2*dd-1)*(-1)**(r-1)*i**(-r)/r
    p3 = sum(sl[l]*T.partial_bell_no1(r, l, T._xargs)/sp.factorial(r) for l in range(1, r//2+1))
    p4 = (2-2*dd)*(-1)**(r-1)*(2*i)**(-r)/r
    s5 = (i-T.aq)**(-r)+(i+T.aq)**(-r)+(i-T.bq)**(-r)+(i+T.bq)**(-r)
    p5 = -(-1)**(r-1)*s5/r
    return sp.radsimp(sp.expand(p12+p3+p4+p5))
t0 = time.time(); n = 0; bad = []
for dd in range(6, 12):
    for key in keys_of_degree(dd):
        Ls = {r: Lnum(r, dd, key[0]) for r in range(1, 9)}
        ser = sp.series(sp.exp(sum(Ls[r]*y**r for r in range(1, 9))), y, 0, 9).removeO()
        e00 = rho(key, 0)
        for delta in range(5, min(8, 2*dd-4)+1):
            val = sp.expand(sp.radsimp(sp.expand(ser.coeff(y, delta))))
            r = rho(key, delta)/e00
            ok = (r.p1 == 0 and r.p3 == 0 and sp.Rational(r.p0.numerator, r.p0.denominator) == sp.re(val) and sp.Rational(r.p2.numerator, r.p2.denominator) == sp.im(val))
            n += 1
            if not ok: bad.append((key, delta))
    print('done d =', dd, f'{time.time()-t0:.0f}s', flush=True)
print(f'rho(5..8) via Bell(L_1..L_8, proved N1..N4) vs recorded data: {n} cells, {len(bad)} bad  ({time.time()-t0:.1f}s)')
print(bad[:8])
