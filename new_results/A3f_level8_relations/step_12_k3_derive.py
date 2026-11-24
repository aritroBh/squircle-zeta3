"""step_12_k3_derive.py: extract K3 equations from the 4-dim symbol solution space;
integrate symmetric tensors to log-cubics G; numeric check + degeneration y->0."""
import sympy as sp
from fractions import Fraction as Fr
from itertools import product as iprod
import mpmath as mp

x, y = sp.symbols('x y')
xi, eta = 1 - x, 1 - y
K4 = [
    (-1, 2, 1, -1, 1), (-1, 1, 2, 1, -1),
    (1, 2, 1, -1, -2), (1, 1, 2, -2, -1),
    (1, 1, 1, 0, 0), (1, 1, 1, -1, -1),
    (-1, 1, 1, -1, 0), (-1, 1, 1, 0, -1),
    (1, 1, 0, 0, 1), (1, 0, 1, 1, 0),
    (-1, 1, 0, -1, 1), (-1, 0, 1, 1, -1),
    (-1, 1, 0, -1, -1), (-1, 0, 1, -1, -1),
    (1, 1, 0, 0, -1), (1, 0, 1, -1, 0),
    (1, 1, 0, 0, 0), (1, 0, 1, 0, 0),
    (-1, 1, 0, -1, 0), (-1, 0, 1, 0, -1),
]
def mon(s, i, j, k, l):
    return sp.Integer(s) * x**i * y**j * xi**k * eta**l
def canon(L):
    L = sp.expand(L)
    s, sn = sp.sstr(L), sp.sstr(-L)
    return min(s, sn)
def factor_vec(expr):
    num, den = sp.fraction(sp.together(expr))
    out = {}
    for e, sgn in ((num, 1), (den, -1)):
        content, facs = sp.factor_list(sp.Poly(e, x, y).as_expr())
        for p_, e_ in facs:
            out[canon(p_)] = out.get(canon(p_), 0) + sgn * e_
    return {k: v for k, v in out.items() if v}
LSET = {}
for s, i, j, k, l in K4:
    m = mon(s, i, j, k, l)
    for d in (factor_vec(m), factor_vec(1 - m)):
        for key in d: LSET[key] = None
LKEYS = sorted(LSET)
LIDX = {k: n for n, k in enumerate(LKEYS)}
NL = len(LKEYS)
def lin_vec(d):
    v = [0] * NL
    for k, e in d.items(): v[LIDX[k]] += e
    return v
syms = []
for s, i, j, k, l in K4:
    m = mon(s, i, j, k, l)
    vm = lin_vec(factor_vec(m)); v1m = lin_vec(factor_vec(1 - m))
    T = {}
    for a_, b_, c_ in iprod(range(NL), repeat=3):
        if v1m[a_] and vm[b_] and vm[c_]:
            T[(a_, b_, c_)] = T.get((a_, b_, c_), 0) + v1m[a_] * vm[b_] * vm[c_]
    syms.append(T)

def combine(cvec):
    S = {}
    for T, c in zip(syms, cvec):
        for k3, v in T.items():
            S[k3] = S.get(k3, 0) + c * v
    return {k3: Fr(v).limit_denominator() for k3, v in S.items() if v}

# solution: last 4 coords free (c16..c19); express all 20 in terms of them
# from step_11: with (c16,c17,c18,c19)=(p,q,r,s):
def sol(p, q, r, s):
    p, q, r, s = Fr(p), Fr(q), Fr(r), Fr(s)
    return [p/6 + q/3 + r/6 - s/6, p/3 + q/6 - r/6 + s/6, p/6 - q/6 + r/6 + s/3,
            -p/6 + q/6 + r/3 + s/6, -p - q, -r - s, -q - r, -p - s,
            -p/2 - q/2 + s/2, -p/2 - q/2 + r/2, -q/2 - r/2 + s/2, -p/2 + r/2 - s/2,
            q/2 - r/2 - s/2, p/2 - r/2 - s/2, -p/2 + q/2 - s/2, p/2 - q/2 - r/2,
            p, q, r, s]

# symmetric (x<->y) subfamily: c16=c17=p, c18=c19=q
for (p, q) in [(1, 0), (0, 1), (6, 6), (6, -6)]:
    cv = sol(p, q, p, q)
    S = combine(cv)
    # check symmetry
    sym_ok = all(S.get((a, b, c), 0) == S.get((b, a, c), 0) == S.get((a, c, b), 0)
                 for (a, b, c) in list(S))
    print(f'(p,q)=({p},{q}): symmetric={sym_ok}, |S|={len(S)}, coeffs={[str(v) for v in cv]}')

# G for a symmetric S: G = (1/6) sum s_abc log|La|log|Lb|log|Lc|
mp.mp.dps = 60
def letter_val(key, xv, yv):
    return abs(float(sp.expand(sp.sympify(key)).subs({x: xv, y: yv})))
def G_val(S, xv, yv):
    tot = mp.mpf(0)
    logs = [mp.log(letter_val(k, xv, yv)) for k in LKEYS]
    for (a, b, c), v in S.items():
        tot += mp.mpf(str(v)) * logs[a] * logs[b] * logs[c] / 6
    return tot
def F_val(cv, xv, yv):
    tot = mp.mpf(0)
    for (s, i, j, k, l), c in zip(K4, cv):
        if c == 0: continue
        m = s * xv**i * yv**j * (1-xv)**k * (1-yv)**l
        tot += mp.mpf(str(c)) * mp.polylog(3, m)
    return tot

for (p, q) in [(1, 0), (0, 1), (6, 6), (6, -6), (1, 1), (1, -1)]:
    cv = sol(p, q, p, q)
    S = combine(cv)
    pt = (mp.mpf('0.2'), mp.mpf('0.3'))
    Fv, Gv = F_val(cv, *pt), G_val(S, *pt)
    print(f'(p,q)=({p},{q}): F-G at (0.2,0.3) = {mp.nstr(Fv-Gv, 12)}  F={mp.nstr(Fv,8)} G={mp.nstr(Gv,8)}')
    # degeneration y -> 0+
    pt2 = (mp.mpf('0.2'), mp.mpf('1e-50'))
    Fv2, Gv2 = F_val(cv, *pt2), G_val(S, *pt2)
    print(f'    degeneration y=1e-50: F={mp.nstr(Fv2,12)} G={mp.nstr(Gv2,12)} diff={mp.nstr(Fv2-Gv2,6)}')
