"""Planner check: the full closed form of M4 = R5 + R''_nice + W3^(4) + T, against M4a's M4_800.txt.

- R5: Codex's closed_form.json (W5), numerical identity.
- R''_nice: M4a w4pp_closed_form.json (= M4b2's proved form).
- W3^(4) = W3' - zeta(3)/8 - pi^2 lambda/24 - lambda^3/12, with W3' from main.tex (proved).
- T = (pi^2 + 6 lambda^2)/24 (proved).
"""
import json, re, sys
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp

sys.argv = [sys.argv[0], sys.argv[1] if len(sys.argv) > 1 else '820']
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_eval_R5 as E          # computes R5 at import (prints its own comparison)

r2, lam, L = E.r2, E.lam, E.L
b = r2 - 1
def fr(s): f = F(s); return mp.mpf(f.numerator) / f.denominator

alias = {'log2': 'lam', 'pi^2': 'pi2', 'Li2(b)': 'Li2(b)', 'Li3(b)': 'Li3(b)', 'Li3(-b)': 'Li3((1-1*sqrt2))',
         'Li3(1/sqrt2)': 'Li3(1/2*sqrt2)'}
W4pp = mp.mpf(0)
for c, name in json.loads((HERE.parent.parent / 'M4a_decomposition' / 'w4pp_closed_form.json').read_text())['closed_form_terms']:
    v = fr(c)
    for f in E.split(name):
        if f == 'sqrt2': v *= r2; continue
        v *= E.atom(alias.get(f, f))
    W4pp += v
W3p = (7 * r2 / 32 - mp.mpf(3) / 16) * mp.zeta(3) + (1 + r2) * mp.polylog(3, b) + (1 - r2) * mp.polylog(3, -b) \
      - 2 * r2 * mp.polylog(3, 1 / r2) + (2 * L - r2 * lam) * mp.polylog(2, b) \
      + (r2 / 12 - mp.mpf(1) / 8) * mp.pi ** 2 * L - (1 + r2) / 24 * mp.pi ** 2 * lam + (1 + r2) / 24 * lam ** 3 \
      + r2 / 4 * lam ** 2 * L - (1 + r2) / 2 * lam * L ** 2 + (5 - r2) / 6 * L ** 3
W34 = W3p - mp.zeta(3) / 8 - mp.pi ** 2 * lam / 24 - lam ** 3 / 12
T = (mp.pi ** 2 + 6 * lam ** 2) / 24
M4 = E.R5 + W4pp + W34 + T
M4_ref = mp.mpf(re.sub(r'\s', '', (HERE.parent.parent / 'M4a_decomposition' / 'M4_800.txt').read_text()).split('=')[-1])
W4pp_ref = mp.mpf(re.sub(r'\s', '', (HERE.parent.parent / 'M4a_decomposition' / 'W4pp_800.txt').read_text()).split('=')[-1])
print("R''_nice - W4''(M4a 800) =", mp.nstr(W4pp - W4pp_ref, 8))
print('M4 closed form =', mp.nstr(M4, 50))
print('M4 closed form - M4_800 (M4a) =', mp.nstr(M4 - M4_ref, 8))
