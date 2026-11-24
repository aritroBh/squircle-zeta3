import sys
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M5a_next_coefficient')
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
from functools import cache
from fractions import Fraction as F
from flint import acb, acb_series, acb_poly, ctx
import mpmath as mp
from qsqrt2 import QS2 as Q

ctx.prec = 300
ctx.cap = 200
N = ctx.cap
mp.mp.dps = 60

def qa(c):
    return acb(str(c.p)) + acb(2).sqrt()*acb(str(c.q))

@cache
def real_series(poles):
    if not poles: return acb_series([1])
    assert poles[-1] != Q(0)
    tail = real_series(poles[1:])
    if poles[0] == Q(0): return acb_series([0]+[tail[j]/j for j in range(1, N)])
    return (tail/acb_series([-qa(poles[0]), 1])).integral()

@cache
def real_g_raw(w, c, e):
    if not w: return acb(1)
    if c == e: return acb(0)
    def taylor(w, c, e):
        if w[-1] == c: return None
        nz = [abs((p-c)/(e-c)) for p in w if p != c]
        if nz and min(nz) <= Q(F(3, 2)): return None
        poles = tuple((p-c)/(e-c) for p in w)
        return acb_poly(real_series(poles).coeffs())(acb(1))
    value = taylor(w, c, e)
    if value is not None: return value
    value = taylor(tuple(reversed(w)), e, c)
    if value is not None: return (-1)**len(w)*value
    mid = (c+e)/2
    return sum((real_g(w[:k], mid, e)*real_g(w[k:], c, mid) for k in range(len(w)+1)), acb(0))

def li_mn(m, n, x, y):
    xf, yf = mp.mpf(float(x)), mp.mpf(float(y))
    return mp.fsum(xf**p * yf**q / (p**m * q**n) for p in range(1, 3000) for q in range(1, p))

for (m, n, x, y) in [(3, 2, F(1,3), F(1,2)), (4, 2, F(1,3), F(1,2)), (3, 3, F(1,4), F(1,3))]:
    xq, yq = Q(x), Q(y)
    w = (Q(0),)*(m-1) + (Q(1)/xq,) + (Q(0),)*(n-1) + (Q(1)/(xq*yq),)
    val = real_g_raw(w, Q(0), Q(1))
    ser = li_mn(m, n, x, y)
    got = mp.mpf(val.real.mid().str(50, radius=False, more=True))
    print(f"Li_{{{m},{n}}}({x},{y}): series={mp.nstr(ser,45)}")
    print(f"   G_0(0^{m-1},1/x,0^{n-1},1/(xy);1) = {mp.nstr(got,45)}, ratio = {mp.nstr(got/ser, 8)}")