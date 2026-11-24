import sys, json
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/M5a_next_coefficient'
sys.path.insert(0, work)
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
from functools import cache
from fractions import Fraction as F
from flint import acb, acb_series, acb_poly, ctx
import mpmath as mp
from qsqrt2 import QS2 as Q, SQ2 as R, A, B

ctx.prec = 400
ctx.cap = 250
N = ctx.cap
mp.mp.dps = 80

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
    return sum((real_g_raw(w[:k], mid, e)*real_g_raw(w[k:], c, mid) for k in range(len(w)+1)), acb(0))

ONE = Q(1); S2I = Q(F(0), F(1, 2))
atoms = [((3,3), S2I, ONE), ((3,3), B, Q(2)), ((4,2), B, S2I)]
vals = json.load(open(f"{work}/n6_values.json"))
keys = ["Li33_s1", "Li33_b2", "Li42_bs"]
for ((m, n), x, y), key in zip(atoms, keys):
    w = (Q(0),)*(m-1) + (ONE/x,) + (Q(0),)*(n-1) + (ONE/(x*y),)
    v = real_g_raw(w, Q(0), ONE)
    got = mp.mpf(v.real.mid().str(70, radius=False, more=True))
    ser = mp.mpf(vals[key])
    print(f"Li_{{{m},{n}}}({x},{y}): |integral - series| = {mp.nstr(abs(got-ser), 4)}  rad={v.real.rad().str(4)}")
with open(f'{work}/progress.log', 'a') as f:
    f.write("step: N6 atom values cross-checked series vs flint iterated-integral evaluator\n")