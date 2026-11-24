# Task (1): independent rewrite of the degree-5 residue polynomials in the binomial basis C(k-6, m-1)
# and exact comparison with e5m_binomial_coeffs.json.  Also exact comparison against the raw kernel
# pickles (kern_k.pkl, k = 6..18) of the whole family.
import json, pickle, sys
from fractions import Fraction
from math import comb
import sympy as sp

BASE = '<project>/research/08_theory/'
RES = BASE + 'TOWN1_tower_large_n/claude_checks/deg5_residue_polynomials.json'
RES_WT6S = BASE + 'WT6S_weight6_symbol_all_k/deg5_residue_polynomials.json'
EJ = BASE + 'WT6S_weight6_symbol_all_k/e5m_binomial_coeffs.json'
KDIR = BASE + 'WTS_weight_k_structure/theory/'

k = sp.symbols('k')
res = json.load(open(RES))
res2 = json.load(open(RES_WT6S))
print('TOWN1 copy == WT6S copy of deg5_residue_polynomials.json :', res == res2)
ej = json.load(open(EJ))

def polyvals(expr):
    """expr: string polynomial in k -> function k -> Fraction (exact)"""
    e = sp.sympify(expr, locals={'k': k})
    P = sp.Poly(e, k)
    assert P.degree() <= 4
    coeffs = [(Fraction(int(c.p), int(c.q))) for c in P.all_coeffs()]  # highest first
    def f(kk):
        kk = Fraction(kk)
        out = Fraction(0)
        for c in coeffs:
            out = out * kk + c
        return out
    return f

# pole labels: pickle pole tuple (r0,r1) = r0 + r1*sqrt2
POLE_TUP = {'0': (0, 0), 'a': (1, 1), '-a': (-1, -1), '-b': (1, -1)}
KEYMAP = {'50': (5, 0), '41': (4, 1), '32': (3, 2)}

# --- (1) binomial rewrite via forward differences: p(x)=sum_m C(x,m-1) Delta^{m-1}p(0), x = k-6
mine = {}   # (key,pole) -> ({'rat':[e1..e5], 'sqrt2':[...]})
for key, poles in res.items():
    for pole, parts in poles.items():
        out = {}
        for part in ('rat', 'sqrt2'):
            if part in parts:
                f = polyvals(parts[part])
            else:
                f = lambda kk: Fraction(0)
            vals = [f(6 + x) for x in range(5)]
            e = []
            for m in range(1, 6):
                n = m - 1
                e.append(sum((-1) ** (n - j) * comb(n, j) * vals[j] for j in range(n + 1)))
            out[part] = e
        mine[(key, pole)] = out

allok = True
nmatch = 0
for (key, pole), out in mine.items():
    ref = ej[f'{key}|{pole}']
    for part in ('rat', 'sqrt2'):
        r = [Fraction(s) for s in ref[part]]
        ok = (r == out[part])
        nmatch += ok
        if not ok:
            allok = False
            print('MISMATCH', key, pole, part, 'mine', [str(x) for x in out[part]], 'ref', ref[part])
print('cells (key,pole,part) compared:', 2 * len(mine), ' exact matches:', nmatch, ' ALL MATCH:', allok)
print('keys in e5m_binomial_coeffs.json:', len(ej), ' my cells:', len(mine))
print('same key set:', set(ej.keys()) == set(f'{a}|{b}' for (a, b) in mine))

# reconstruction check: sum_m C(k-6,m-1) e_m == polynomial value, for many k (exact)
bad = 0
for (key, pole), out in mine.items():
    for part in ('rat', 'sqrt2'):
        if part in res[key][pole]:
            f = polyvals(res[key][pole][part])
        else:
            f = lambda kk: Fraction(0)
        for kk in list(range(6, 53)) + [100, 1000]:
            val = sum(comb(kk - 6, m - 1) * out[part][m - 1] for m in range(1, 6))
            if val != f(kk):
                bad += 1
print('reconstruction mismatches over k in 6..52,100,1000 (all 15 cells x 2 parts):', bad)

# E_{5,5} coefficient = 24 * (k^4 coefficient)
ok24 = True
for (key, pole), out in mine.items():
    for part in ('rat', 'sqrt2'):
        if part in res[key][pole]:
            P = sp.Poly(sp.sympify(res[key][pole][part], locals={'k': k}), k)
            c4 = P.coeff_monomial(k ** 4)
            c4 = Fraction(int(c4.p), int(c4.q))
        else:
            c4 = Fraction(0)
        if out[part][4] != 24 * c4:
            ok24 = False
print('e_{ij,r,5} == 24 * k^4-coefficient, all cells/parts:', ok24)

# --- cross-check against the raw kernels: coefficient of 1/(v-r) in key (i,j) at weight-k kernel
import os
tot = 0
badk = []
for kk in range(6, 19):
    path = KDIR + f'kern_{kk}.pkl'
    if not os.path.exists(path):
        continue
    K = pickle.load(open(path, 'rb'))
    for key, kt in KEYMAP.items():
        for pole, pt in POLE_TUP.items():
            real = K[kt]['real']
            # pickle coefficient = (c0,c1,c2,c3): c0 + c1 sqrt2 + i (c2 + c3 sqrt2)
            c = real.get((Fraction(pt[0]), Fraction(pt[1])))
            assert c is not None, (kk, kt, pole)
            c0, c1, c2, c3 = [Fraction(x) for x in c]
            rat = sum(comb(kk - 6, m - 1) * mine[(key, pole)]['rat'][m - 1] for m in range(1, 6))
            sq = sum(comb(kk - 6, m - 1) * mine[(key, pole)]['sqrt2'][m - 1] for m in range(1, 6))
            tot += 1
            if (rat, sq, c2, c3) != (c0, c1, 0, 0):
                badk.append((kk, key, pole))
    # no other poles / no imaginary parts / no ipos, ineg, const in these keys
    for key, kt in KEYMAP.items():
        d = K[kt]
        extra = set(d['real'].keys()) - set((Fraction(p[0]), Fraction(p[1])) for p in POLE_TUP.values())
        if extra or d['ipos'] or d['ineg'] or any(Fraction(x) != 0 for x in d['const']):
            print('NOTE extra data at k', kk, key, extra, d['ipos'], d['ineg'], d['const'])
print('pickle cells compared (k=6..18):', tot, ' mismatches vs binomial-basis e (rat,sqrt2,0,0):', len(badk), badk[:5])

# dump my own e's (independently derived) for use in later steps
json.dump({f'{key}|{pole}': {p: [str(x) for x in v] for p, v in out.items()} for (key, pole), out in mine.items()},
          open('my_e5m.json', 'w'), indent=1)
print('wrote my_e5m.json')
