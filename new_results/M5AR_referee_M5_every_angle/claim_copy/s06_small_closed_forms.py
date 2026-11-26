"""s06: the small parts of M_5 written with A = artanh c, L = log(1-c^2) and the dilogarithmic function
   D(c) = H_{0,1}(c) + H_{0,-1}(c) - H_{-1,1}(c) + H_{1,-1}(c) = Li2(c) - Li2(-c) - H_{-1,1}(c) + H_{1,-1}(c)
(c^e He[0,1] - c^e He[-1,1] at odd e).  Exact check: Taylor coefficients c^-6..c^200 of the closed forms vs the
term lists (exact rationals).  Both sides are Laurent series convergent on 0<|c|<1, and the closed forms below are
just a regrouping of the term lists, so the check to c^200 is a check of the regrouping (each regrouping is
term-by-term: He[1] = 2A at odd e, = -L at even e)."""
import json, re
from fractions import Fraction as Fr
import g8_allN as G
NMAX = 200
def series_of_terms(fname):
    out = {}
    for nm, v in json.load(open(fname)):
        q = Fr(v)
        if nm == '1': out[0] = out.get(0, 0) + q; continue
        m = re.match(r'c\^(-?\d+)$', nm)
        if m: out[int(m.group(1))] = out.get(int(m.group(1)), 0) + q; continue
        m = re.match(r'c\^(-?\d+) He\[(.*)\]', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(2).split(','))
        co = G.hpl_coeffs(w, NMAX + 6)
        for k in range(1, NMAX + 6):
            if (k + e) % 2 == 0 and co[k]:
                out[k + e] = out.get(k + e, 0) + 2*q*co[k]
    return {k: v for k, v in out.items() if v and k <= NMAX}
def series_of_terms_plain(fname):
    out = {}
    for nm, v in json.load(open(fname)):
        q = Fr(v)
        if nm == '1': out[0] = out.get(0, 0) + q; continue
        m = re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(2).split(','))
        for k, x in enumerate(G.hpl_coeffs(w, NMAX + 6)):
            if x and k + e <= NMAX: out[k + e] = out.get(k + e, 0) + q*x
    return {k: v for k, v in out.items() if v}
A = {k: Fr(1, k) for k in range(1, NMAX + 7, 2)}                 # artanh c
L = {k: Fr(-2, k) for k in range(2, NMAX + 7, 2)}                # log(1-c^2) = -sum c^{2m}/m
Dw = {}
for w, s in (((0, 1), 1), ((0, -1), 1), ((-1, 1), -1), ((1, -1), 1)):
    for k, v in enumerate(G.hpl_coeffs(w, NMAX + 6)):
        if v: Dw[k] = Dw.get(k, 0) + s*v
def lin(*pairs):
    out = {}
    for laur, ser in pairs:     # laur: {e: q} Laurent polynomial; ser: series dict (None = 1)
        for e, q in laur.items():
            if ser is None: out[e] = out.get(e, 0) + q; continue
            for k, v in ser.items(): out[k + e] = out.get(k + e, 0) + q*v
    return {k: v for k, v in out.items() if v and k <= NMAX}
F = Fr
forms = {
 # [z5]M5 = -5/8 + 3/(8c^2) - 3(1-c^2)^2/(8c^3) A
 'z5': lin(({0: F(-5, 8), -2: F(3, 8)}, None), ({-3: F(-3, 8), -1: F(6, 8), 1: F(-3, 8)}, A)),
 # [z2 z3]M5 = -1/8 + 3/(8c^2) - (1-c^2)(3+c^2)/(8c^3) A
 'z2*z3': lin(({0: F(-1, 8), -2: F(3, 8)}, None), ({-3: F(-3, 8), -1: F(2, 8), 1: F(1, 8)}, A)),
 # G = [z2^2]M5 = 7/16 - 3/(4c^2) - (16 - 19c^2 + c^4)/(32c^4) L + (7 - 5c^2)/(16c^3) A + (c^2+3)(c^2-1)/(32c^3) D
 'z2^2': lin(({0: F(7, 16), -2: F(-3, 4)}, None), ({-4: F(-16, 32), -2: F(19, 32), 0: F(-1, 32)}, L),
             ({-3: F(7, 16), -1: F(-5, 16)}, A), ({-3: F(-3, 32), -1: F(2, 32), 1: F(1, 32)}, Dw)),
}
for nm, f in forms.items():
    print(f"[{nm}]M5 closed form == term list, Laurent coefficients c^-6..c^{NMAX}, exact:", f == series_of_terms(f"hpl_fit_k5_{nm}.json"))
# Psi = [z2]M3 (g8_terms_k3_partz2.json) = 1/2 + (1/(2c^2) - 3/4) L - A/(2c) + (1-c^2)/(4c) D
Psi_form = lin(({0: F(1, 2)}, None), ({-2: F(1, 2), 0: F(-3, 4)}, L), ({-1: F(-1, 2)}, A), ({-1: F(1, 4), 1: F(-1, 4)}, Dw))
print("Psi closed form == g8_terms_k3_partz2.json (plain terms), exact to c^200:", Psi_form == series_of_terms_plain("g8_terms_k3_partz2.json"))
# the pi^4 part of M5 with real zeta values: z4 [z4]M5 + z2^2 [z2^2]M5 = pi^4 ([z4]/90 + [z2^2]/36) = pi^4 (Psi/90 + 7G/180)
import sympy as sp
c = sp.Symbol('c')
def lp(d): return sp.factor(sum(sp.Rational(v.numerator, v.denominator)*c**e for e, v in d.items()))
cf = {'1': ({0: F(1, 2)}, {0: F(7, 16), -2: F(-3, 4)}), 'L': ({-2: F(1, 2), 0: F(-3, 4)}, {-4: F(-1, 2), -2: F(19, 32), 0: F(-1, 32)}),
      'A': ({-1: F(-1, 2)}, {-3: F(7, 16), -1: F(-5, 16)}), 'D': ({-1: F(1, 4), 1: F(-1, 4)}, {-3: F(-3, 32), -1: F(1, 16), 1: F(1, 32)})}
tot = {}
for key, (p, g) in cf.items():
    d = {}
    for e, v in p.items(): d[e] = d.get(e, 0) + v/90
    for e, v in g.items(): d[e] = d.get(e, 0) + 7*v/180
    tot[key] = d
    print(f"   pi^4 part, coefficient of {key}: pi^4 * ({lp(d)})")
ser = lin((tot['1'], None), (tot['L'], L), (tot['A'], A), (tot['D'], Dw))
z4 = series_of_terms("hpl_fit_k5_z4.json"); g = series_of_terms("hpl_fit_k5_z2^2.json")
chk = {k: z4.get(k, 0)/90 + g.get(k, 0)/36 for k in set(z4) | set(g)}
print("pi^4 part as displayed == [z4]M5/90 + [z2^2]M5/36 (term lists), exact to c^200:", ser == {k: v for k, v in chk.items() if v})
print("pi^4 part elementary (no D)?", not tot['D'] or all(v == 0 for v in tot['D'].values()))
