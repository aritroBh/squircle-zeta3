# s11b: L_parsed.pkl vs the PROVED construction tow2lib.Lr (from R_1..R_4 [prove_S12/S34] + algebraic pieces); the sqrt2 constant
# s5 is simplified once per r; the d,m-dependence is then compared by exact rational-function cancellation.
import sys, time
sys.dont_write_bytecode = True
from collections import Counter
import sympy as sp
import tow2lib as T2, tow3lib as T3
d, m = T3.d, T3.m; I = sp.I
t0 = time.time()
Lp = T3._load_L()
def Lr_fast(r):
    p12 = (2*d-1) * (-1)**(r-1) * I**(-r) / r
    p3 = sp.Integer(0)
    for l in range(1, r//2+1):
        p3 += T2.sl[l] * T2.partial_bell_no1(r, l, T2._xargs) / sp.factorial(r)
    p4 = (2-2*d) * (-1)**(r-1) * (2*I)**(-r) / r
    s5 = ((I-T2.aq)**(-r) + (I+T2.aq)**(-r) + (I-T2.bq)**(-r) + (I+T2.bq)**(-r))
    s5c = sp.nsimplify(sp.simplify(sp.radsimp(sp.expand(s5))))
    assert s5c.is_number and not s5c.has(sp.sqrt(2)), s5c
    p5 = -(-1)**(r-1) * s5c / r
    return p12 + p3 + p4 + p5
Lf = {}
for r in range(1, 7):
    Lf[r] = Lr_fast(r)
    diff = sp.cancel(sp.together(Lf[r] - Lp[r]))
    print(f'  L({r}): cached L_parsed - proved construction == 0 : {diff == 0}  ({round(time.time()-t0,1)}s)', flush=True)
def parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from parts(rem-p, p, cur+[p])
for dl in (5, 6):
    tot = sp.Integer(0)
    for part in parts(dl, 1, []):
        C = Counter(part); coef = sp.Integer(1)
        for sz, mult in C.items(): coef *= (Lf[sz]**mult)/sp.factorial(mult)
        tot += coef
    diff = sp.cancel(sp.together(tot - T3.arho(dl)))
    print(f'  a_{dl}: Bell(proved-construction L_1..L_{dl}) - tow3lib.arho({dl}) == 0 : {diff == 0}  ({round(time.time()-t0,1)}s)', flush=True)
