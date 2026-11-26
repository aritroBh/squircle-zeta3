"""s00: test of the fit machinery m5fit.py on k=4 (known answer: Psi_2 = g8_terms_k4_partz2.json, R_4 = ..._part1.json)."""
import sys, json, time
sys.setrecursionlimit(20000)
import sympy as sp
import g8_allN as G, m5fit as F
from fractions import Fraction as Fr
n = sp.Symbol('n')
P = G.mu_side(4)
def conv(E): return {w: (sp.Poly(c.numer.as_expr(), n, domain='QQ'), sp.Poly(c.denom.as_expr(), n, domain='QQ')) for w, c in E.items()}
for zm, W, fn in [((1,0,0), 3, 'g8_terms_k4_partz2.json'), ((0,0,0), 5, 'g8_terms_k4_part1.json')]:
    t0 = time.time()
    nunk = 2 + 5*sum(3**(L-1) for L in range(1, W+1))
    sol, cols, info = F.fit(conv(P[zm]), W, range(-3, 2), nunk + 60, nprimes=2)
    ref = {}
    for nm, v in json.load(open(fn)):
        if nm == '1': ref[('const', 0)] = Fr(v); continue
        e = int(nm.split()[0][2:]); w = tuple(int(t) for t in nm.split('[')[1].rstrip(']').split(','))
        if w[-1] == -1:   # normalise representative: c^e(H_w + s H_-w) = s * c^e(H_-w + s H_w)
            s = 1 if (e + sum(1 for a in w if a)) % 2 == 0 else -1
            w = tuple(-a for a in w); v = Fr(v)*s
        ref[('He', e, w)] = ref.get(('He', e, w), 0) + Fr(v)
    print(G.zmon_name(zm), "fit == G8 data file:", sol == {k: v for k, v in ref.items() if v}, f"({len(sol)} terms, {time.time()-t0:.0f}s)")
