cwd = "<project>/research/08_theory/G2_M4alpha"
code = '''"""step_10.py -- assemble B4 and validate.

B4(c) = R4(c) + z2*Psi2(c) - z3*Psi(c) + z4*F1(c) + z2^2*F2(c), with
  R4   = hpl_fit_k4_(0,0,0).json  (109 terms, weights 2..5; step_09)
  Psi2 = hpl_fit_k4_(1,0,0).json  (22 terms, weights 1..3; step_08)
  Psi  = G1 hpl_fit_k3_(1,0,0).json (B3's z2 part; step_08)
  F1,F2 = b4_pi4_part.json (elementary artanh combos; step_07)

Checks:
 (1) EXACT: full combination's Taylor series == exact target series (all five zeta parts) through
     c^1280, including the vanishing of all odd coefficients.
 (2) NUMERICAL: B4 evaluated from its own HPL Taylor series (order 7000, far beyond the fitted
     c^1280) vs corner_series M4 values at six angles (values4_*.txt, ~90 digits), and vs the
     800-digit squircle value at c = 1/sqrt(2).
"""
import sys, json, re, time
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp
import g2_fit as gf

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
ORDER = 1280
_c = sp.symbols('c')


def load_terms(files):
    out = []
    for name in files:
        for nm, v in json.load(open(name)):
            if nm == '1':
                out.append(('1', (), 0, Fr(v))); continue
            if nm == 'c':
                out.append(('c', (), 1, Fr(v))); continue
            m = re.match(r'c\\^(-?\\d+) (He?)\\[(.*)\\]', nm)
            e = int(m.group(1))
            w = tuple(int(t) for t in m.group(3).split(',')) if m.group(3) else ()
            out.append(('even' if m.group(2) == 'He' else 'hpl', w, e, Fr(v)))
    return out


def laurent_terms(cof):
    """rational function of c (Laurent polynomial) -> list of (Fraction q, int e) with q c^e."""
    num, den = sp.fraction(sp.together(sp.sympify(cof)))
    num = sp.Poly(sp.expand(num), _c); den = sp.Poly(sp.expand(den), _c)
    assert len(den.terms()) == 1, f'non-monomial denominator {den}'
    (dn,), dc = den.terms()[0]
    out = []
    for (nn,), nc in num.terms():
        out.append((Fr(sp.Rational(nc))/Fr(sp.Rational(dc)), nn - dn))
    return out


def pi4_terms():
    """F1, F2 as HPL term lists: R(c)*artanh(c) + S(c); artanh = (H_1 + H_-1)/2."""
    data = json.load(open('b4_pi4_part.json'))
    out = {}
    for fname in ('F1', 'F2'):
        ts = []
        for cof, kind in data[fname]:
            for q, e in laurent_terms(cof):
                if kind == 'artanh':
                    ts.append(('hpl', (1,), e, q/2))
                    ts.append(('hpl', (-1,), e, q/2))
                else:
                    ts.append(('1c', (), e, q))
        out[fname] = ts
    return out


def term_series(kind, w, e, order):
    if kind == '1':
        return [Fr(1)] + [Fr(0)]*order
    if kind == 'c':
        return [Fr(0), Fr(1)] + [Fr(0)]*(order - 1)
    if kind == '1c':
        return [Fr(1) if j == e else Fr(0) for j in range(order + 1)]
    if kind == 'hpl':
        s = gf.hpl_series(tuple(w), order)
        return [s[j - e] if 0 <= j - e <= order else Fr(0) for j in range(order + 1)]
    if kind == 'even':
        return gf.even_series(tuple(w), e, order)
    raise ValueError(kind)


def combo_series(terms, order):
    T = [Fr(0)]*(order + 1)
    for kind, w, e, q in terms:
        s = term_series(kind, w, e, order)
        for j in range(order + 1):
            if s[j]:
                T[j] += q*s[j]
    return T


def exact_check():
    tab = json.load(open('flat_exact_N640_K4_formula.json'))['4']
    p4 = pi4_terms()
    parts = {
        '(0, 0, 0)': load_terms(['hpl_fit_k4_(0,0,0).json']),
        '(1, 0, 0)': load_terms(['hpl_fit_k4_(1,0,0).json']),
        '(0, 1, 0)': [(k, w, e, -q) for k, w, e, q in load_terms([G1 + 'hpl_fit_k3_(1,0,0).json'])],
        '(0, 0, 1)': p4['F1'],
        '(2, 0, 0)': p4['F2'],
    }
    ok_all = True
    for key, ts in parts.items():
        T = combo_series(ts, ORDER)
        ok = True
        for nn in range(1, 641):
            if T[2*nn] != Fr(tab[str(nn)].get(key, '0')):
                ok = False; print(f'  part {key}: mismatch at n={nn}'); break
        odd_ok = all(T[m] == 0 for m in range(1, ORDER + 1, 2))
        print(f'  part {key}: {len(ts)} terms; even match through c^1280: {ok}; odd vanish: {odd_ok}')
        ok_all = ok_all and ok and odd_ok
    print('(1) EXACT through c^1280, all five parts:', ok_all)
    return ok_all


# ---------------- numerical ----------------
NUMORD = 7000
_num_cache = {}


def num_hpl_series(w):
    w = tuple(w)
    if w in _num_cache:
        return _num_cache[w]
    N = NUMORD + 3
    if len(w) == 1:
        a = w[0]
        s = [mp.mpf(0)]*(N + 1)
        for m in range(1, N + 1):
            s[m] = mp.mpf(1)/m if a == 1 else -mp.mpf((-1)**m)/m
        _num_cache[w] = s
        return s
    inner = num_hpl_series(w[1:])
    a = w[0]
    out = [mp.mpf(0)]*(N + 1)
    if a == 0:
        for m in range(1, N + 1):
            out[m] = inner[m]/m
    else:
        acc = mp.mpf(0)
        for m in range(0, N):
            acc = acc*a + inner[m]
            out[m + 1] = acc/(m + 1)
    _num_cache[w] = out
    return out


def num_even_series(w, e):
    s1 = num_hpl_series(w)
    w2 = tuple(-a for a in w)
    s2 = num_hpl_series(w2) if w2 != w else s1
    sgn = 1 if (e + sum(1 for a in w if a)) % 2 == 0 else -1
    return [s1[m] + sgn*s2[m] for m in range(NUMORD + 1)]


def eval_terms(terms, c):
    tot = mp.mpf(0)
    for kind, w, e, q in terms:
        qq = mp.mpf(q.numerator)/mp.mpf(q.denominator)
        if kind == '1':
            tot += qq
            continue
        if kind == 'c':
            tot += qq*c
            continue
        if kind == '1c':
            tot += qq*c**e
            continue
        s = num_hpl_series(w) if kind == 'hpl' else num_even_series(w, e)
        v = mp.mpf(0)
        cp = c**max(0, -e)
        for m in range(max(0, -e), NUMORD + 1):
            v += s[m]*cp
            cp *= c
        tot += qq*v*(c**e if e >= 0 else 1)
    return tot


def b4_value(c, verbose=False):
    R4 = load_terms(['hpl_fit_k4_(0,0,0).json'])
    Psi2 = load_terms(['hpl_fit_k4_(1,0,0).json'])
    Psi = [(k, w, e, -q) for k, w, e, q in load_terms([G1 + 'hpl_fit_k3_(1,0,0).json'])]
    z2 = mp.pi**2/6; z3 = mp.zeta(3)
    p4 = json.load(open('b4_pi4_part.json'))
    v = mp.mpf(0)
    for fname, mul in (('F1', mp.pi**4/90), ('F2', mp.pi**4/36)):
        for cof, kind in p4[fname]:
            cv = mp.mpf(0)
            for q, e in laurent_terms(cof):
                cv += (mp.mpf(q.numerator)/q.denominator)*c**e
            v += mul*cv*(mp.atanh(c) if kind == 'artanh' else 1)
    vR = eval_terms(R4, c); v2 = eval_terms(Psi2, c); v3 = eval_terms(Psi, c)
    if verbose:
        print(f'    R4={mp.nstr(vR, 30)}\\n    z2*Psi2={mp.nstr(z2*v2, 30)}\\n    -z3*Psi={mp.nstr(z3*v3, 30)}\\n    pi4={mp.nstr(v, 30)}')
    return vR + z2*v2 + z3*v3 + v


if __name__ == '__main__':
    t0 = time.time()
    ok = exact_check()
    print(f'(1) exact check done ({time.time()-t0:.0f}s)', flush=True)
    mp.mp.dps = 120
    ANG = {'pi_6': mp.pi/6, 'pi_4': mp.pi/4, 'pi_3': mp.pi/3, 'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'gen1p3': mp.mpf(13)/10}
    M4SQ = mp.mpf(open('<project>/research/08_theory/M4a_decomposition/M4_800.txt').read().strip())
    for tag, a in ANG.items():
        c = mp.cos(a/2)
        val = b4_value(c)
        ref = mp.mpf(open(f'values4_{tag}.txt').read().split('M4 =')[1].strip())
        d = abs(val - ref)
        extra = ''
        if tag == 'pi_2':
            extra = f'; |B4 - M4_800| = {mp.nstr(abs(val - M4SQ), 3)}'
        print(f'(2) {tag:>7} (c={mp.nstr(c, 6)}): |B4 - M4 quadrature| = {mp.nstr(d, 3)}{extra}', flush=True)
'''
with open(cwd + "/step_10.py", "w") as f:
    f.write(code)
import ast; ast.parse(code)
print("written")
