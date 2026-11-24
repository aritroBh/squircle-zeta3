code = r'''"""step_11_k3_symbol.py: derive Kummer's two-variable trilogarithm equation (K3) via the symbol method.

Ansatz: same 20 monomial arguments m_i(x,y) as Kummer's K4 (functions of x,y,xi=1-x,eta=1-y).
Find c in Q^20, c != 0, such that  S = sum_i c_i * (1-m_i) x m_i x m_i   (tensor symbol,
letters = non-constant irreducible factors; constants like -1 and rationals dropped)
is a totally SYMMETRIC tensor. Then S = symbol of an explicit log-cubic G(x,y), and
  sum_i c_i Li3(m_i) = G   on D0 (modulo constants fixed by degeneration y->0+).
"""
import sympy as sp
from fractions import Fraction as Fr
from itertools import product as iprod

x, y = sp.symbols('x y')
xi, eta = 1 - x, 1 - y

# the 20 Kummer arguments as (sign, i, j, k, l): m = sign * x^i y^j xi^k eta^l
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

letters = [x, y, xi, eta]
def canon(L):
    """canonical letter key: primitive integer polynomial up to sign"""
    L = sp.expand(L)
    cont, fac = sp.factor_list(sp.Poly(L, x, y).as_expr()) if True else (1, L)
    s = sp.sstr(L)
    sn = sp.sstr(-L)
    return min(s, sn)

def factor_vec(expr):
    """expr (rational in x,y) -> dict letter_key -> exponent (non-constant factors only)."""
    num, den = sp.fraction(sp.together(expr))
    out = {}
    for e, sgn in ((num, 1), (den, -1)):
        content, facs = sp.factor_list(sp.Poly(e, x, y).as_expr())
        for p_, e_ in facs:
            out[canon(p_)] = out.get(canon(p_), 0) + sgn * e_
    return {k: v for k, v in out.items() if v}

# build letter list
LSET = {}
for s, i, j, k, l in K4:
    m = mon(s, i, j, k, l)
    for d in (factor_vec(m), factor_vec(1 - m)):
        for key in d:
            LSET[key] = None
LKEYS = sorted(LSET)
LIDX = {k: n for n, k in enumerate(LKEYS)}
NL = len(LKEYS)
print('letters:', NL)
for k in LKEYS: print('   ', k)

def lin_vec(d):
    v = [0] * NL
    for k, e in d.items():
        v[LIDX[k]] = v[LIDX[k]] + e
    return v

# symbol of each argument: (1-m) x m x m as sparse tensor dict (a,b,c) -> coef
syms = []
for s, i, j, k, l in K4:
    m = mon(s, i, j, k, l)
    vm = lin_vec(factor_vec(m))
    v1m = lin_vec(factor_vec(1 - m))
    T = {}
    for a_, b_, c_ in iprod(range(NL), repeat=3):
        if v1m[a_] and vm[b_] and vm[c_]:
            T[(a_, b_, c_)] = T.get((a_, b_, c_), 0) + v1m[a_] * vm[b_] * vm[c_]
    syms.append(T)

# unknowns c_1..c_20; equations: sum c_i T_i symmetric under S3 (swap 1<->2 and 2<->3)
cs = sp.symbols('c0:20')
eqs = {}
for T, c in zip(syms, cs):
    for (a, b, d), v in T.items():
        eqs[(a, b, d)] = eqs.get((a, b, d), 0) + c * v
conds = set()
for (a, b, d), e in eqs.items():
    conds.add(sp.expand(e - eqs.get((b, a, d), 0)))
    conds.add(sp.expand(e - eqs.get((a, d, b), 0)))
conds = [sp.expand(e_) for e_ in conds if e_ != 0]
print('symmetry conditions:', len(conds), flush=True)
sol = sp.linsolve(conds, cs)
print('solution set:', sol)
'''
open('step_11_k3_symbol.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_11_k3_symbol.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout); print('STDERR:', r.stderr[-3000:])
