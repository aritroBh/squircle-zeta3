
import os, json, datetime
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
step = r'''"""WT6S step A: E_{5,m} exactly as Q(sqrt2)-combinations of N_{(i,j),r} from deg5 JSON.

res_{ij,r}(k) = rat(k) + sqrt2*sq(k) (quartic, k=6..52). Rewrite in the binomial basis
C(k-6,m-1), m=1..5: res(k) = sum_m C(k-6,m-1) e_{ij,r,m}. Then
  W6(k) = sum_{ij,r} res_{ij,r}(k) N_{(i,j),r} = sum_m C(k-6,m-1) E_{5,m},
  E_{5,m} = sum_{ij,r} e_{ij,r,m} N_{(i,j),r}.
Checks: (A1) binomial rewrite exact at k=6,7,10,21 vs the polynomials; (A2) the leading
coefficient combination equals the LEAD-law constant c(6) numerically (Task A quadrature)."""
import json
from fractions import Fraction as F
import sympy as sp

k = sp.Symbol('k')
KEYS = ["50", "41", "32"]
POLES = ["0", "a", "-a", "-b"]

def parse_poly(s):
    return sp.expand(sp.sympify(s)) if s is not None else sp.Integer(0)

def binomial_coeffs(poly):
    """Exact coefficients c[m], m=0..4, with poly(k) = sum_m c[m]*C(k-6,m) (c in Q)."""
    poly = sp.Poly(sp.expand(poly), k)
    assert poly.degree() <= 4
    # binomial basis matrix: C(k-6,m) = sum_j M[j,m] k^j  -> solve
    out = []
    rem = sp.expand(poly)
    for m in range(4, -1, -1):
        binom = sp.expand(sp.binomial(k - 6, m))
        lc = sp.LC(sp.Poly(binom, k))
        cdeg = sp.Poly(rem, k).degree() if rem != 0 else -1
        if cdeg < m:
            out.append(sp.Integer(0)); continue
        assert cdeg == m
        coeff = sp.Rational(sp.Poly(rem, k).nth(m)) / lc
        out.append(coeff)
        rem = sp.expand(rem - coeff * binom)
    assert rem == 0
    out.reverse()   # index m
    return out

def main():
    dd = json.load(open('deg5_residue_polynomials.json'))
    # e[(key, pole)] = {'rat': [e0..e4], 'sqrt2': [e0..e4]}  (rational lists)
    e = {}
    for key in KEYS:
        for pole in POLES:
            cell = dd[key].get(pole, {})
            rat = parse_poly(cell.get('rat'))
            sq = parse_poly(cell.get('sqrt2'))
            e[(key, pole)] = {'rat': [str(c) for c in binomial_coeffs(rat)],
                              'sqrt2': [str(c) for c in binomial_coeffs(sq)]}
    json.dump({f"{a}|{b}": e[(a, b)] for a in KEYS for b in POLES},
              open('e5m_binomial_coeffs.json', 'w'), indent=1)
    # A1: binomial rewrite reproduces the polynomials exactly at k=6,7,10,21
    ok = True
    for key in KEYS:
        for pole in POLES:
            cell = dd[key].get(pole, {})
            for part in ['rat', 'sqrt2']:
                poly = parse_poly(cell.get(part))
                for kv in [6, 7, 10, 21, 52]:
                    val = sum(sp.Rational(e[(key, pole)][part][m]) * sp.binomial(kv - 6, m)
                              for m in range(5))
                    ok &= (sp.expand(val - poly.subs(k, kv)) == 0)
    print("A1: binomial rewrite reproduces residue polys at k=6,7,10,21,52:", ok, flush=True)
    assert ok
    # E_{5,m} coefficient table (rational + sqrt2 parts), printed
    for m in range(1, 6):
        print(f"--- E_{{5,{m}}} = sum e_{{ij,r,{m}}} N_{{(i,j),r}} ---")
        for key in KEYS:
            i, j = int(key[0]), int(key[1])
            for pole in POLES:
                pr = sp.Rational(e[(key, pole)]['rat'][m - 1])
                qs = sp.Rational(e[(key, pole)]['sqrt2'][m - 1])
                if pr or qs:
                    cstr = f"({pr}" + (f" + {qs}*sqrt2)" if qs else ")")
                    print(f"   {cstr} * N_{{({i},{j}),{pole}}}")
    # A2 (exact part): leading coefficient of each residue = e_5/24; print sum for c(6)
    lead = {}
    for key in KEYS:
        for pole in POLES:
            cell = dd[key].get(pole, {})
            lr = sp.LC(sp.Poly(parse_poly(cell.get('rat')), k)) if parse_poly(cell.get('rat')) != 0 else sp.Integer(0)
            lq = sp.LC(sp.Poly(parse_poly(cell.get('sqrt2')), k)) if parse_poly(cell.get('sqrt2')) != 0 else sp.Integer(0)
            lead[(key, pole)] = (str(lr), str(lq))
            # consistency: e5 = 24 * leading coeff (since C(k-6,4) ~ k^4/24)
            assert sp.Rational(e[(key, pole)]['rat'][4]) == 24 * lr
            assert sp.Rational(e[(key, pole)]['sqrt2'][4]) == 24 * lq
    json.dump({f"{a}|{b}": lead[(a, b)] for a in KEYS for b in POLES},
              open('e5m_lead_coeffs.json', 'w'), indent=1)
    print("A2 exact: e_5 = 24 * leading-coefficient for every cell: PASS", flush=True)
    print("c(6) = sum over (ij,r) lead-coeff * N_{(i,j),r}; numeric check in step_E.", flush=True)
    with open('progress.log', 'a') as f:
        f.write('step_A: deg5 residues rewritten in binomial basis C(k-6,m-1); A1 exact at k=6,7,10,21,52; e5=24*lead PASS\n')
    print("saved e5m_binomial_coeffs.json, e5m_lead_coeffs.json")

if __name__ == '__main__':
    main()
'''
open(os.path.join(WD, 'step_03.py'), 'w').write(step)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'step_03.py'], capture_output=True, text=True, cwd=WD)
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
