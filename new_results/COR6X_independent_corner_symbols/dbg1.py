
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, json, re
path = "<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"
data = json.load(open(path))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def symbol_Omega(c0):
    tot = {}
    for w, coeff in TERMS:
        tot = sym_add(tot, sym_scale(sym_H(w, c0), coeff))
    return tot
def sym_Li4(z):
    return sym_scale(sym_G((0,0,0,1), z), -1)

S_om = sym_expanded(symbol_Omega(sp.Rational(1,3)))
S_rhs = sym_expanded(sym_add(
    sym_scale(sym_Li4(-8), sp.Rational(-11,81)),
    sym_scale(sym_Li4(-3), sp.Rational(-2)),
    sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(20,3))))
print("=== S(Omega(1/3)) expanded, primes {2,3}: (monomials over (2,3))")
for k, v in sorted(S_om.items(), key=str):
    print("   ", k, v)
print("=== S(RHS) expanded:")
for k, v in sorted(S_rhs.items(), key=str):
    print("   ", k, v)
print("=== diff:")
diff = sym_sub(S_om, S_rhs)
for k, v in sorted(diff.items(), key=str):
    print("   ", k, v)
print("delta22(diff) =", delta22(diff))
print("delta13-primitive check: ker of both cobrackets?")
