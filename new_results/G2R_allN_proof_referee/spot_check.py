"""spot_check.py -- G2R referee, task 4: independent spot check, n = 1..40, exact in Q, per zeta monomial.

mu_{4,n}: from the DEFINITION, via G1's flat_exact.py output table (flat_run/flat_exact_N40_K4.json,
produced by running their unmodified script; integral-representation evaluation, Bell/cumulant exact).
[c^N] of B4: built here with an OWN parser and an OWN exact HPL Taylor evaluator (nothing imported
from step_10/step_11/g2_fit/hpl/zsum). He[w] means c^e (H_w + (-1)^{e+nz(w)} H_{-w}).
Also cross-checks: mu4_formula.pkl (the proof's mu input, step_05's symbolic formula) evaluated at
harmonic numbers == flat_exact for all n <= 40; and the odd coefficients and low orders c^-3..c^1.
"""
import os, json, re, pickle, time
from fractions import Fraction as Fr
from functools import lru_cache
import sympy as sp

BASE = os.path.dirname(os.path.abspath(__file__))
G2 = "<project>/research/08_theory/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
t0 = time.time()
NMAX = 40
ORD = 2*NMAX          # highest even power probed

# ---------- own exact HPL Taylor evaluator ----------
@lru_cache(maxsize=None)
def hc(w, m):
    """[c^m] H_w(c), letters {0,1,-1}, last letter nonzero; exact Fraction."""
    w = tuple(w)
    if m <= 0:
        return Fr(0)
    if len(w) == 1:
        a = w[0]
        return Fr(a**(m - 1), m)                      # H_1: 1/m ; H_-1: (-1)^(m+1)/m
    a, rest = w[0], w[1:]
    if a == 0:
        return hc(rest, m)/m
    return sum(((a**(m - 1 - i))*hc(rest, i) for i in range(1, m)), Fr(0))/m

def nz(w):
    return sum(1 for a in w if a != 0)

def neg(w):
    return tuple(-a for a in w)

# ---------- own parsers ----------
def load_hpl_fit(path):
    """-> list of (kind, w, e, q): kind in {'1','c','He','H'}."""
    out = []
    for nm, v in json.load(open(path)):
        q = Fr(v)
        if nm == '1':
            out.append(('1', (), 0, q)); continue
        if nm == 'c':
            out.append(('c', (), 1, q)); continue
        m = re.match(r'c\^(-?\d+) (He|H)\[(.*)\]', nm)
        out.append((m.group(2), tuple(int(t) for t in m.group(3).split(',')), int(m.group(1)), q))
    return out

def load_F(fname):
    """b4_pi4_part.json: (7c/16 - 5/(8c) + 3/(16c^3)) artanh + 9/16 - 3/(16c^2); artanh = (H_1+H_-1)/2."""
    data = json.load(open(os.path.join(G2, 'b4_pi4_part.json')))[fname]
    c = sp.symbols('c')
    out = []
    for cof, kind in data:
        num, den = sp.fraction(sp.together(sp.sympify(cof)))
        num = sp.Poly(sp.expand(num), c); den = sp.Poly(sp.expand(den), c)
        assert len(den.terms()) == 1
        (dn,), dc = den.terms()[0]
        for (ne,), nc in num.terms():
            q = Fr(sp.Rational(nc))/Fr(sp.Rational(dc)); e = ne - dn
            if kind == 'artanh':
                out.append(('H', (1,), e, q/2))
                out.append(('H', (-1,), e, q/2))
            else:
                out.append(('1c', (), e, q))
    return out

def coeff_of(terms, N):
    """[c^N] of the combination, exact."""
    tot = Fr(0)
    for kind, w, e, q in terms:
        if kind == '1':
            if N == 0: tot += q
        elif kind == 'c':
            if N == 1: tot += q
        elif kind == '1c':
            if N == e: tot += q
        elif kind == 'H':
            m = N - e
            if m >= 1: tot += q*hc(w, m)
        elif kind == 'He':
            m = N - e
            if m >= 1:
                sgn = 1 if (e + nz(w)) % 2 == 0 else -1
                tot += q*(hc(w, m) + sgn*hc(neg(w), m))
    return tot

PARTS = {
    '(0, 0, 0)': load_hpl_fit(G2 + 'hpl_fit_k4_(0,0,0).json'),
    '(1, 0, 0)': load_hpl_fit(G2 + 'hpl_fit_k4_(1,0,0).json'),
    '(0, 1, 0)': [(k, w, e, -q) for k, w, e, q in load_hpl_fit(G1 + 'hpl_fit_k3_(1,0,0).json')],
    '(0, 0, 1)': load_F('F1'),
    '(2, 0, 0)': load_F('F2'),
}
ref = json.load(open(os.path.join(BASE, 'flat_run/flat_exact_N40_K4.json')))['4']

ok_even = ok_odd = ok_low = True
resid = []
for part, terms in PARTS.items():
    for nv in range(1, NMAX + 1):
        got = coeff_of(terms, 2*nv)
        want = Fr(ref[str(nv)].get(part, '0'))
        if got != want:
            ok_even = False; resid.append((part, nv, got, want))
    for m in range(1, ORD + 1, 2):
        if coeff_of(terms, m) != 0:
            ok_odd = False; resid.append((part, f'odd c^{m}', coeff_of(terms, m), 0))
    for N in range(-3, 2):
        if coeff_of(terms, N) != 0:
            ok_low = False; resid.append((part, f'low c^{N}', coeff_of(terms, N), 0))
print(f"(4a) [c^(2n)] B4 part == flat_exact mu_{{4,n}} part, exact, n=1..{NMAX}, all five parts: {ok_even}")
print(f"(4b) odd coefficients c^1..c^{ORD-1} vanish, all parts: {ok_odd}")
print(f"(4c) low orders c^-3..c^1 vanish, all parts: {ok_low}")
if resid:
    print("residuals:", resid[:6])

# ---------- cross-check the proof's mu input: mu4_formula.pkl vs the definition ----------
terms = pickle.load(open(G2 + 'mu4_formula.pkl', 'rb'))
n = sp.Symbol('n')
KEY = {(0,0,0): '(0, 0, 0)', (1,0,0): '(1, 0, 0)', (0,1,0): '(0, 1, 0)', (0,0,1): '(0, 0, 1)', (2,0,0): '(2, 0, 0)'}
def H(r, m):
    return sum((Fr(1, i**r) for i in range(1, m + 1)), Fr(0))
okF = True
for nv in range(1, NMAX + 1):
    hv = [H(r, nv - 1) for r in range(1, 5)] + [H(r, 2*nv - 1) for r in range(1, 5)]
    acc = {k: Fr(0) for k in KEY.values()}
    for mon, cf in terms:
        val = Fr(sp.Rational(sp.sympify(cf).subs(n, nv)))
        for x, e in zip(hv, mon[:8]):
            if e:
                val *= x**e
        acc[KEY[(mon[8], mon[9], mon[10])]] += val
    for k in KEY.values():
        if acc[k] != Fr(ref[str(nv)].get(k, '0')):
            okF = False; print("   mu4_formula mismatch", nv, k, acc[k], ref[str(nv)].get(k, '0'))
print(f"(4d) mu4_formula.pkl (proof's mu input) == flat_exact definition, all n<= {NMAX}, all parts: {okF}")
print(f"SPOT CHECK DONE ({time.time()-t0:.0f}s):", ok_even and ok_odd and ok_low and okF)
