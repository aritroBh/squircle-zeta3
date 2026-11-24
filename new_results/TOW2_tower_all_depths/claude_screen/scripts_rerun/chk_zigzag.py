# Independent data test of the zigzag top law e(2s,s)/e00 = (-1)^s A_{s+1}/Pi_s, e(2s+1,s) = -i(d-s-2) e(2s,s)
# off-diagonal keys only; Newton coefficient via finite differences of the recorded towers (k = d..d+s <= 14)
import sympy as sp
from towlib import *
from math import comb
xx = sp.symbols('xx')
ser = sp.series(sp.sech(xx)+sp.tanh(xx), xx, 0, 12).removeO().expand()
A = {n: abs(sp.Rational(ser.coeff(xx, n))*sp.factorial(n)) for n in range(12)}
print('zigzag A_n:', [A[n] for n in range(1, 9)])
def e_from_data(key, delta, s):
    i, j = key; d = i+j
    base = d
    return sum((u(i, j, base+t, delta).scale(Fr((-1)**(s-t)*comb(s, t))) for t in range(0, s+1)), ZERO)
tested = {}; bad = []
for s in range(1, 7):
    for d in range(s+2, 15):
        if d + s > 14: continue
        for key in keys_of_degree(d):
            i, j = key
            if i == j: continue
            e00 = e_from_data(key, 0, 0)
            Pi = 1
            for jj in range(1, s+1): Pi *= (2*d-2*jj-1)
            expect = sp.Rational((-1)**s * int(A[s+1]), Pi)
            ee = e_from_data(key, 2*s, s)/e00
            ok = (ee.p1 == 0 and ee.p3 == 0 and ee.p2 == 0 and sp.Rational(ee.p0.numerator, ee.p0.denominator) == expect)
            tested.setdefault(s, []).append((d, key, ok))
            if not ok: bad.append((s, key, ee))
            # odd top
            if 2*s+1 <= 2*d-4:
                eo = e_from_data(key, 2*s+1, s)/e00
                exp_o = -sp.I*(d-s-2)*expect
                ok2 = (eo.p1 == 0 and eo.p3 == 0 and sp.Rational(eo.p0.numerator, eo.p0.denominator) == sp.re(exp_o) and sp.Rational(eo.p2.numerator, eo.p2.denominator) == sp.im(exp_o))
                if not ok2: bad.append((s, key, 'odd', eo))
for s, L in tested.items():
    print(f's={s}: even-top tested for degrees', sorted(set(x[0] for x in L)), 'cells', len(L), 'all ok:', all(x[2] for x in L))
print('bad:', bad[:6])
