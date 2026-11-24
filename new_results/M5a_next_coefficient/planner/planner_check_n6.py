"""Planner: check M5a's N6 against W6's (3,3)/(4,2) pairings with my own pairing code (planner_w6_pairing.py helpers).
Li_{m,n}(x,y) = G(0^{m-1}, 1/x, 0^{n-1}, 1/(xy); 1), symbol via w4_symbol.g_symbol(ZERO, word, ONE)."""
import sys
from pathlib import Path
src = Path(__file__).with_name('planner_w6_pairing.py').read_text()
src = src[:src.index('# weight-5 control')]
exec(compile(src, 'helpers', 'exec'))
def li_mn(m, n, x, y):
    w = (ZERO,) * (m - 1) + (ONE / x,) + (ZERO,) * (n - 1) + (ONE / (x * y),)
    return g_symbol(ZERO, w, ONE)
b = R - 1
N6 = add(scale(li_mn(3, 3, b, ONE), Q(-12) + 209 * R), scale(li_mn(3, 3, R / 2, ONE), Q(-160) + 80 * R),
         scale(li_mn(3, 3, b, Q(2)), Q(F(329, 4)) * R), scale(li_mn(4, 2, b, R / 2), 120 * R))
S6 = symbol_of(load(TH / 'M5a_next_coefficient' / 'gpl_w6.tsv'))
print('N6 pairings :', [fmt(v) for v in tests6(N6)])
print('W6 pairings :', [fmt(v) for v in tests6(S6)])
print('W6 - N6     :', [fmt(v) for v in tests6(add(S6, scale(N6, -1)))])
