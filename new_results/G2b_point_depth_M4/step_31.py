
import sys
from fractions import Fraction as Fr, Fraction as F
BASE = "<project>/research/08_theory"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")

# --- COR6 side: symbol of J((0,0,0,0,2); 1) over basis {2, 1+sqrt2} (exponent vectors)
from cor6_sym import sym_J
def vec_q(el):  # el rational -> exponent vector in basis {2, 1+sqrt2}; rationals are powers of 2 only here
    el = abs(el)
    assert el.denominator & (el.denominator - 1) == 0 and el.numerator & (el.numerator - 1) == 0
    import math
    e = int(math.log2(el.numerator)) - int(math.log2(el.denominator))
    return (Fr(e), Fr(0))
S_cor6 = sym_J((0,0,0,0,2), 1, vec_q, a0=0)
print("COR6 sym_J((0,0,0,0,2);1):", {w: str(c) for w, c in S_cor6.items()})

# --- M4c side: g_symbol(0, (0,0,0,0,2), 1)
from w4_symbol import g_symbol, ZERO, ONE, li_symbol
S_m4c = g_symbol(ZERO, (ZERO,ZERO,ZERO,ZERO, ONE*2), ONE)
print("M4c g_symbol(G_0((0,0,0,0,2);1)):", {w: (str(c.p), str(c.q)) for w, c in S_m4c.items()})
print("M4c -li_symbol(5, 1/2):", {w: (str(c.p), str(c.q)) for w, c in li_symbol(5, F(1,2)).items()})
