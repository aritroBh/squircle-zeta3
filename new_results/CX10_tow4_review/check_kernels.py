"""Independent full partial-fraction and WT6S checks, using exact Fractions."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import pickle
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
KEYS = ((5, 0), (4, 1), (3, 2))
K = {k: pickle.load(open(ROOT/'WTS_weight_k_structure'/'theory'/f'kern_{k}.pkl', 'rb'))
     for k in range(6, 12)}
WT = ROOT/'WT6S_weight6_symbol_all_k'
TABLE = json.loads((WT/'e5m_binomial_coeffs.json').read_text())
POLYS = json.loads((WT/'deg5_residue_polynomials.json').read_text())
POLES = {'0': (F(0), F(0)), 'a': (F(1), F(1)), '-a': (F(-1), F(-1)),
         'b': (F(-1), F(1)), '-b': (F(1), F(-1))}
ZERO = (F(0),)*4

def flat(pf):
    out = {('real', pole, c): F(v) for pole, cs in pf['real'].items() for c, v in enumerate(cs)}
    out.update({(part, n, c): F(v) for part in ('ipos', 'ineg')
                for n, cs in enumerate(pf[part]) for c, v in enumerate(cs)})
    out.update({('const', 0, c): F(v) for c, v in enumerate(pf['const'])})
    return out

def combine(pfs, weights):
    terms = list(map(flat, pfs))
    return {n: sum(w*t.get(n, 0) for w, t in zip(weights, terms))
            for n in set().union(*terms)}

nreal = 0
kvar = sp.Symbol('k')
for key in KEYS:
    pfs = [K[k][key] for k in range(6, 12)]
    assert not any(combine(pfs, [1, -5, 10, -10, 5, -1]).values()), key
    nreal += 4*len(set().union(*(pf['real'] for pf in pfs)))
    for k, pf in zip(range(6, 12), pfs):
        assert all(F(v) == 0 for part in ('ipos', 'ineg') for cs in pf[part] for v in cs), (key, k)
        assert all(F(v) == 0 for v in pf['const']), (key, k)
    print('PASS', key, 'complete K6 = 5K7-10K8+10K9-5K10+K11; +/-i and constants zero k=6..11', flush=True)
    tag = ''.join(map(str, key))
    assert set(POLYS[tag]) == {'0', 'a', '-a', '-b'}
    for order in range(5):
        weights = [(-1)**(order-j)*comb(order, j) for j in range(order+1)]
        target = {'real': {}, 'ipos': [], 'ineg': [], 'const': ZERO}
        for name in POLYS[tag]:
            cell = TABLE[f'{tag}|{name}']
            target['real'][POLES[name]] = (F(cell['rat'][order]), F(cell['sqrt2'][order]), F(0), F(0))
        assert not any(combine(pfs[:order+1]+[target], weights+[-1]).values()), (key, order)
        print('PASS', key, f'Delta^{order} K6 == WT6S E_5,{order+1} coefficient kernel', flush=True)
    for name, parts in POLYS[tag].items():
        cell = TABLE[f'{tag}|{name}']
        for part in ('rat', 'sqrt2'):
            # Product form avoids unevaluated symbolic binomial at negative k.
            newton = sum(sp.Rational(cell[part][j]) * sp.prod(kvar-6-h for h in range(j))/sp.factorial(j)
                         for j in range(5))
            assert sp.expand(newton-sp.sympify(parts.get(part, '0'))) == 0, (key, name, part)
    for k, pf in zip(range(6, 12), pfs):
        expected = {POLES[name]: tuple(F(sp.sympify(parts.get(part, '0')).subs(kvar, k))
                    for part in ('rat', 'sqrt2'))+(F(0), F(0)) for name, parts in POLYS[tag].items()}
        target = {'real': expected, 'ipos': [], 'ineg': [], 'const': ZERO}
        assert not any(combine([pf, target], [1, -1]).values()), (key, k)
    print('PASS', key, 'WT6S power polynomials == Newton table symbolically and all kernels k=6..11', flush=True)
print(f'PASS {nreal}/48 real QSI components, 18 full PF zero checks, 15 difference kernels, 24 polynomial identities', flush=True)
