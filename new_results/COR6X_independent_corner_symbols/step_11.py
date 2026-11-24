
# step_11.py -- task (2): symbol-level localization of the bug at c0=1/2
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, json, re
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
c0 = sp.Rational(1, 2)
S_om = {}
for w, cf in TERMS:
    S_om = sym_add(S_om, sym_scale(sym_H(w, c0), cf))
S_om = sym_expanded(S_om)

def sym_Li4(z): return sym_scale(sym_G((0,0,0,1), z), -1)
S_five = sym_expanded(sym_add(
    sym_scale(sym_Li4(-2), sp.Rational(-21)),
    sym_scale(sym_Li4(sp.Rational(1,3)), sp.Rational(-65,8)),
    sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(13,2)),
    sym_scale(sym_Li4(sp.Rational(1,2)), sp.Rational(137,6)),
    sym_scale(sym_Li4(sp.Rational(1,9)), sp.Rational(39,64))))

# symbols of the product terms P
S_log2 = {(sp.Integer(2),): 1}; S_log3 = {(sp.Integer(3),): 1}
S_Li3_13 = sym_expanded(sym_scale(sym_G((0,0,1), sp.Rational(1,3)), -1))   # Li_3(1/3) = -G(0,0,1;1/3)
S_Li3_m12 = sym_expanded(sym_scale(sym_G((0,0,1), sp.Rational(-1,2)), -1)) # Li_3(-1/2)
S_Li2_13 = sym_expanded(sym_scale(sym_G((0,1), sp.Rational(1,3)), -1))     # Li_2(1/3)
def sh(*args):
    r = args[0]
    for a in args[1:]:
        r = sym_shuffle(r, a)
    return r
S_P = {}
S_P = sym_add(S_P, sym_scale(sh(S_Li3_13, S_log2), sp.Rational(15,2)))
S_P = sym_add(S_P, sym_scale(sh(S_Li3_m12, S_log2), sp.Rational(37,2)))
S_P = sym_add(S_P, sym_scale(sh(S_Li3_m12, S_log3), sp.Rational(-5,4)))
S_P = sym_add(S_P, sh(S_Li2_13, S_Li2_13))
S_P = sym_add(S_P, sym_scale(sh(S_Li2_13, S_log2, S_log2), sp.Rational(-23,4)))
S_P = sym_add(S_P, sym_scale(sh(S_Li2_13, S_log2, S_log3), sp.Rational(7,4)))
S_P = sym_add(S_P, sh(S_Li2_13, S_log3, S_log3))
S_P = sym_add(S_P, sym_scale(sh(S_log2, S_log2, S_log2, S_log2), sp.Rational(-211,144)))
S_P = sym_add(S_P, sym_scale(sh(S_log2, S_log2, S_log2, S_log3), sp.Rational(101,24)))
S_P = sym_add(S_P, sym_scale(sh(S_log2, S_log2, S_log3, S_log3), sp.Rational(-15,8)))
S_P = sym_add(S_P, sym_scale(sh(S_log2, S_log3, S_log3, S_log3), sp.Rational(-35,24)))
S_P = sym_add(S_P, sym_scale(sh(S_log3, S_log3, S_log3, S_log3), sp.Rational(113,192)))
S_P = sym_expanded(S_P)

T = sym_sub(S_om, sym_add(S_five, S_P))
print("S(Omega(1/2)) - S(fiveLi4) - S(P): terms =", len(T))
for k, v in sorted(T.items(), key=str):
    print("   ", k, v)
print()
print("S(Omega) terms:", len(S_om), " S(five) terms:", len(S_five), " S(P) terms:", len(S_P))
