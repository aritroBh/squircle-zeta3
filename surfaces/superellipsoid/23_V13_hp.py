# 23_V13_hp.py -- V1 or V3' = V3 - V2 alone at high precision (for the LLL attempt in 22_lll_vertex_try.py).
import mpmath as mp, sys, json
which = sys.argv[1]; DPS = int(sys.argv[2])
mp.mp.dps = DPS + 25
r2 = mp.sqrt(2); L = mp.log(1+r2)
def qu(u): return mp.sqrt(u*u + (1-u)**2)
def I(f): return mp.quad(f, [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1])
def safe_old(f):
    def g(u):
        u = mp.mpf(u)
        with mp.workdps(mp.mp.dps + (int(-mp.log10(u)) + 10 if u < mp.mpf('0.01') else 0)):
            v = f(u)
        return +v
    return g
def safe(f):
    # no extra precision: near u=0 the tanh-sinh weights are O(u), so the O(eps/u) rounding error of the
    # 0/0 brackets contributes only O(eps) to the integral
    return f
if which == "V1":
    val = I(safe(lambda u: qu(u)*mp.log1p(3*u*u-2*u)*mp.log(u)/u))
else:   # V3' = int [2(1-u) atanh((1-u)/(sqrt2 q)) - 2L] log u du/u
    val = I(safe(lambda u: (2*(1-u)*mp.atanh((1-u)/(r2*qu(u))) - 2*L)*mp.log(u)/u))
json.dump({which: mp.nstr(val, DPS)}, open("%s_%d.json" % (which, DPS), "w"))
print(which, mp.nstr(val, 50))
