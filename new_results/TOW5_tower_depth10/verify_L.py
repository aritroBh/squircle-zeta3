"""verify_L.py -- L cache == tow2lib-style L(r) construction, symbolic, r = 1..10.
For r = 9, 10 the construction uses sl_5 (tow5lib, PROVED from R_5). sqrt2-aware reduction:
numerator of the difference, sqrt2 -> w, even/odd w-parity parts both == 0.
"""
import sys, os, time, pickle
sys.dont_write_bytecode = True
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work)
import sympy as sp
import tow2lib as T2
import tow5lib as T5   # sets T2.sl[5]
d, m = T2.d, T2.m
w = sp.Symbol('w')
def qsqrt2_is_zero(expr):
    g = sp.together(expr)
    num = sp.fraction(g)[0]
    num = sp.expand(num.subs({sp.sqrt(2): w}))
    P = sp.Poly(num, w)
    a = sp.Integer(0); b = sp.Integer(0)
    for (e,), c in P.terms():
        if e % 2 == 0: a += c * 2**(e//2)
        else: b += c * 2**((e-1)//2)
    return sp.expand(a) == 0 and sp.expand(b) == 0
ok = True
t0 = time.time()
for r in range(1, 11):
    t1 = time.time()
    diff = T5._L[r] - T2.L(r)
    z = qsqrt2_is_zero(diff)
    ok &= z
    print(f'r={r}: L cache == tow2lib-style L(r) in Q(i,sqrt2)(d,m): {z}  ({round(time.time()-t1,1)}s)', flush=True)
print('VERDICT L (r=1..10):', ok, f'total {round(time.time()-t0,1)}s', flush=True)
assert ok
