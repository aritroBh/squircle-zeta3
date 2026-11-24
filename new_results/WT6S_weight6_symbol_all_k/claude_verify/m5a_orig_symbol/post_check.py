"""Cross-checks on the outputs of run_cases.py (exact Q(sqrt2) arithmetic with the original QS2 class)."""
import json
import sys
from fractions import Fraction as F
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'orig'))
from w4_symbol import Q  # noqa: E402

WT = ROOT.parent.parent
KEY = ['50', '41', '32']
POLES = ['0', 'a', '-a', '-b']


def load(name):
    raw = json.loads((ROOT / name).read_text())
    return {k: [Q(F(c['p']), F(c['q'])) for c in v] for k, v in raw.items()}


def eq(u, v):
    return all(a == b for a, b in zip(u, v))


def show(vec):
    return '(' + ', '.join(str(x) for x in vec) + ')'


ok_all = True
e5 = json.loads((WT / 'e5m_binomial_coeffs.json').read_text())

# ---- 1. recompute e_{ij,r,m} from deg5_residue_polynomials.json by forward differences (independent of step_03)
import sympy as sp
k = sp.Symbol('k')
rp = json.loads((WT / 'deg5_residue_polynomials.json').read_text())
bad = 0
ncell = 0
for key, poles in rp.items():
    for pole, parts in poles.items():
        for part, tag in (('rat', 'rat'), ('sqrt2', 'sqrt2')):
            poly = sp.Poly(sp.sympify(parts[part], locals={'k': k}), k) if part in parts else sp.Poly(0, k)
            vals = [poly.eval(6 + j) for j in range(5)]
            for m in range(1, 6):
                em = sum((-1) ** (m - 1 - j) * comb(m - 1, j) * vals[j] for j in range(m))   # Delta^{m-1} at k=6
                want = sp.Rational(e5[f'{key}|{pole}'][tag][m - 1])
                if em != want:
                    bad += 1
                    print('MISMATCH e', key, pole, tag, m, em, want)
            ncell += 1
print(f'[1] forward-difference recomputation of all e_(ij,r,m) from deg5_residue_polynomials.json: cells*parts={ncell}, mismatches={bad}')
ok_all &= bad == 0
# poles present in e5m_binomial_coeffs but not in polynomials (or vice versa)
print('    pole/key sets equal:', {f'{k_}|{p}' for k_, ps in rp.items() for p in ps} == set(e5.keys()))

# ---- 2. linearity: sum of 12 single-cell pairing vectors == pairing of E_{5,m} fed as a kernel
N = load('out_group_N.json')
E = load('out_group_E.json')
for m in range(1, 6):
    tot = [Q(0)] * 4
    for ck, v in e5.items():
        key, pole = ck.split('|')
        c = Q(F(v['rat'][m - 1]), F(v['sqrt2'][m - 1]))
        pv = N[f'N{key}|{pole}']
        tot = [t + c * x for t, x in zip(tot, pv)]
    same = eq(tot, E[f'E5,{m}'])
    ok_all &= same
    print(f'[2] E5,{m}: direct-kernel {show(E[f"E5,{m}"])}  | linear combo of 12 N-cells {show(tot)}  equal={same}')

# ---- 3. compare with WT6S's published table
tab = json.loads((WT / 'e5m_pairings.json').read_text())['pairings_pq']
for m in range(1, 6):
    ref = [Q(F(p), F(q)) for p, q in tab[f'E5,{m}']]
    same = eq(ref, E[f'E5,{m}'])
    ok_all &= same
    print(f'[3] E5,{m} vs WT6S e5m_pairings.json: equal={same}')
nt = json.loads((WT / 'n6_symbols_pairings.json').read_text())
print('    n6_symbols_pairings.json keys sample:', list(nt.keys())[:4] if isinstance(nt, dict) else type(nt))

# ---- 4. kernels at k=6..9 built directly from the residue polynomials vs binomial combination of pv(E_m)
for kk in (6, 7, 8, 9):
    tot = [Q(0)] * 4
    for m in range(1, 6):
        b = comb(kk - 6, m - 1)
        tot = [t + b * x for t, x in zip(tot, E[f'E5,{m}'])]
    K = load('out_group_K.json')[f'kpoly{kk}']
    same = eq(tot, K)
    ok_all &= same
    print(f'[4] k={kk}: direct kernel {show(K)}  | sum_m C(k-6,m-1) pv(E5,m) {show(tot)}  equal={same}')

# ---- 5. W6(k) pairing polynomials claimed in result.md section C
sq = lambda p, q: Q(F(p), F(q))
for kk in (6, 7, 8, 9):
    kf = F(kk)
    claim = [sq(0, -3 * kf * kf + F(57, 2) * kf - 17), sq(-5 * kf * kf + F(95, 2) * kf - F(285, 2), 0),
             sq(0, 15), sq(6 * kf * kf - 57 * kf + 126, 0)]
    K = load('out_group_K.json')[f'kpoly{kk}']
    same = eq(claim, K)
    ok_all &= same
    print(f'[5] k={kk}: result.md polynomial P(k) {show(claim)} equals direct kernel pairing: {same}')
print('ALL POST CHECKS PASS' if ok_all else 'SOME POST CHECK FAILED')
