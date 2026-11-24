
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
step = r'''"""WT6S step C: pairings of W6(k) as polynomials in k; verify vs direct kernels.

For k in {6,7,10,21,52} build the W6 symbol from the deg5 residue polynomials evaluated
at k, compute the pairing vector, and compare EXACTLY with
  P(k) = sum_m C(k-6,m-1) * pv(E_{5,m}).
Also print the pairing polynomials in the power basis in k (exact)."""
import json
from fractions import Fraction as F
import sympy as sp
import qsqrt2
from w4_symbol import Q, R, ban_float, self_check
import wt6s_core as C6
from wt6s_core import (ZEROQ, symbol_of_kernel_terms, pairing_vector6, add, log)

k = sp.Symbol('k')
KEYS5 = ["50", "41", "32"]
POLES5 = ["0", "a", "-a", "-b"]

def parse_poly(s):
    return sp.expand(sp.sympify(s)) if s is not None else sp.Integer(0)

def omega_terms_at(kv):
    dd = json.load(open('deg5_residue_polynomials.json'))
    out = {}
    for key in KEYS5:
        i, j = int(key[0]), int(key[1])
        terms = []
        for pole in POLES5:
            cell = dd[key].get(pole, {})
            rat = parse_poly(cell.get('rat')).subs(k, kv)
            sq = parse_poly(cell.get('sqrt2')).subs(k, kv)
            if rat or sq:
                terms.append([pole, [str(sp.Rational(rat)), str(sp.Rational(sq))]])
        out[f"{i},{j}"] = terms
    return out

def w6_pairings_at(kv):
    om = omega_terms_at(kv)
    S = {}
    for ij in [(5, 0), (4, 1), (3, 2)]:
        S = add(S, symbol_of_kernel_terms(om[f"{ij[0]},{ij[1]}"], *ij))
    return pairing_vector6(S)

def parse_q(s):
    s = s.replace(' ', '')
    if '*sqrt2' in s:
        a, bq = s.split('*sqrt2')
        return Q(F(a or '0'), F(bq or '0'))
    return Q(F(s), F(0))

def main():
    ban_float(); self_check()
    ePv = json.load(open('e5m_pairings.json'))
    epv = {m: [parse_q(v) for v in ePv[f"E5,{m}"]] for m in range(1, 6)}
    # verify at k = 6,7,10,21,52
    ok = True
    for kv in [6, 7, 10, 21, 52]:
        direct = w6_pairings_at(kv)
        poly = [sum(sp.binomial(kv - 6, m - 1) * epv[m][t] for m in range(1, 6))
                for t in range(4)]
        match = all(direct[t] == poly[t] for t in range(4))
        ok &= match
        log(f"k={kv}: direct={[str(v) for v in direct]}")
        log(f"      poly  ={[str(v) for v in poly]}  match={match}")
        assert match, kv
    log(f"C: pairing polynomial reproduces direct kernel pairings at k=6,7,10,21,52: {ok}")
    # pairing polynomials in the power basis
    labels = ['33', '42a', '42b', '42c']
    print()
    for t in range(4):
        # P_t(k) = sum_m C(k-6,m-1) epv[m][t], as Q(sqrt2) poly in k
        expr = sum(sp.binomial(k - 6, m - 1) * (epv[m][t].p + sp.sqrt(2) * epv[m][t].q)
                   for m in range(1, 6))
        expr = sp.expand(expr)
        rat = sp.expand(expr).subs(sp.sqrt(2), 0)
        sqr = sp.Rational(str(sp.expand(expr).coeff(sp.sqrt(2), 1)))
        log(f"pairing[{labels[t]}](k) = ({rat}) + ({sqr})*sqrt2")
    json.dump({'note': 'P_t(k)=sum_m C(k-6,m-1) pv(E5m)[t]; verified k=6,7,10,21,52',
               'e5m_pairings': ePv}, open('w6k_pairing_polynomials.json', 'w'), indent=1)
    with open('progress.log', 'a') as f:
        f.write('step_C: W6(k) pairing polynomial == direct kernel pairings at k=6,7,10,21,52 EXACT; power-basis polys printed\n')
    log("saved w6k_pairing_polynomials.json")

if __name__ == '__main__':
    main()
'''
open(os.path.join(WD, 'step_05.py'), 'w').write(step)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'step_05.py'], capture_output=True, text=True, cwd=WD)
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
