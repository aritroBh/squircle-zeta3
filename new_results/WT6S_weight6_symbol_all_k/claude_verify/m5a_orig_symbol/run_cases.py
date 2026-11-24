"""Feed kernel residue vectors through M5a's ORIGINAL weight-six symbol code (orig/w6_symbol.py, verbatim copy).

Usage: run_cases.py GROUP   with GROUP in {E, K, N}
  E : E_{5,m}, m=1..5, from WT6S e5m_binomial_coeffs.json (residue e_{ij,r,m} at pole r on key (i,j))
  K : kernels at k = 5,6,7,8,9 built DIRECTLY from deg5_residue_polynomials.json (independent of the binomial rewrite)
  N : the twelve single-cell N_{(i,j),r} (coefficient 1) -> linearity cross-check
For each case: omega5.json is written in cases/<name>/, the original w6_symbol.HERE is pointed at that folder and the
ORIGINAL w6_expression() / symbol() / pair33() / pair42() are called unchanged.  symbol6.json (same format as M5a's) is
written so that the verbatim planner6_pairing.py can be run on it afterwards.
"""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'orig'))
import w6_symbol as W                      # noqa: E402  (ORIGINAL code, verbatim copy)
from w4_symbol import ban_float, self_check, Q, R  # noqa: E402

WT = ROOT.parent.parent                    # WT6S folder
KEY = {'50': '5,0', '41': '4,1', '32': '3,2'}
POLES = ['0', 'a', '-a', '-b']
ban_float()
self_check()


def omega_from_cells(cells):
    """cells: {(key, pole): (p, q)} -> omega dict in M5a's omega5.json format {"i,j": [[pole, [p, q]], ...]}"""
    om = {'3,2': [], '4,1': [], '5,0': []}
    for (key, pole), (p, q) in cells.items():
        om[KEY[key]].append([pole, [str(F(p)), str(F(q))]])
    return om


def run_case(name, cells, write_symbol=True):
    d = ROOT / 'cases' / name
    d.mkdir(parents=True, exist_ok=True)
    om = omega_from_cells(cells)
    (d / 'omega5.json').write_text(json.dumps(om, indent=1) + '\n')
    W.HERE = d                              # original w6_expression() reads HERE/'omega5.json'
    expr = W.w6_expression()
    S = W.symbol(expr)
    p33 = W.pair33(S)
    p42 = W.pair42(S)
    if write_symbol:
        data = {'W6_pair33': str(p33), 'W6_pair42': [str(v) for v in p42], 'alphabet': ['log2', 'log(1+sqrt2)'],
                'symbol': {''.join(map(str, w)): [str(c.p), str(c.q)] for w, c in sorted(S.items())}}
        (d / 'symbol6.json').write_text(json.dumps(data, indent=1) + '\n')
    print(f"{name:10s} gpl_terms={len(expr):5d} sym_words={len(S):5d}  (3,3)={p33}  (4,2)=({p42[0]}, {p42[1]}, {p42[2]})",
          flush=True)
    return p33, p42


def qs(v):
    return {'p': str(v.p), 'q': str(v.q)}


if __name__ == '__main__':
    grp = sys.argv[1]
    out = {}
    if grp == 'E':
        e = json.loads((WT / 'e5m_binomial_coeffs.json').read_text())
        for m in range(1, 6):
            cells = {}
            for ck, v in e.items():
                key, pole = ck.split('|')
                cells[(key, pole)] = (F(v['rat'][m - 1]), F(v['sqrt2'][m - 1]))
            p33, p42 = run_case(f'E5{m}', cells)
            out[f'E5,{m}'] = [qs(p33)] + [qs(v) for v in p42]
    elif grp == 'K':
        import sympy as sp
        rp = json.loads((WT / 'deg5_residue_polynomials.json').read_text())
        k = sp.Symbol('k')
        for kk in (5, 6, 7, 8, 9):
            cells = {}
            for key, poles in rp.items():
                for pole, parts in poles.items():
                    r = sp.sympify(parts['rat'], locals={'k': k}).subs(k, kk)
                    s = sp.sympify(parts['sqrt2'], locals={'k': k}).subs(k, kk) if 'sqrt2' in parts else sp.Integer(0)
                    cells[(key, pole)] = (F(int(sp.numer(r)), int(sp.denom(r))), F(int(sp.numer(s)), int(sp.denom(s))))
            p33, p42 = run_case(f'kpoly{kk}', cells)
            out[f'kpoly{kk}'] = [qs(p33)] + [qs(v) for v in p42]
    elif grp == 'N':
        for key in ('50', '41', '32'):
            for pole in POLES:
                p33, p42 = run_case(f'N{key}_{pole}', {(key, pole): (F(1), F(0))}, write_symbol=False)
                out[f'N{key}|{pole}'] = [qs(p33)] + [qs(v) for v in p42]
    (ROOT / f'out_group_{grp}.json').write_text(json.dumps(out, indent=1) + '\n')
    print('saved', f'out_group_{grp}.json', flush=True)
