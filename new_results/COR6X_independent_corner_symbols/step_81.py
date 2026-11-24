t = r'''
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
S_log2 = {(sp.Integer(2),): 1}; S_log3 = {(sp.Integer(3),): 1}
S_Li3_13 = sym_expanded(sym_scale(sym_G((0,0,1), sp.Rational(1,3)), -1))
S_Li3_m12 = sym_expanded(sym_scale(sym_G((0,0,1), sp.Rational(-1,2)), -1))
S_Li2_13 = sym_expanded(sym_scale(sym_G((0,1), sp.Rational(1,3)), -1))
def sh(*args):
    r = args[0]
    for a in args[1:]: r = sym_shuffle(r, a)
    return r
S_P = sym_add(
    sym_scale(sh(S_Li3_13, S_log2), sp.Rational(15,2)),
    sym_scale(sh(S_Li3_m12, S_log2), sp.Rational(37,2)),
    sym_scale(sh(S_Li3_m12, S_log3), sp.Rational(-5,4)),
    sh(S_Li2_13, S_Li2_13),
    sym_scale(sh(S_Li2_13, S_log2, S_log2), sp.Rational(-23,4)),
    sym_scale(sh(S_Li2_13, S_log2, S_log3), sp.Rational(7,4)),
    sh(S_Li2_13, S_log3, S_log3),
    sym_scale(sh(S_log2, S_log2, S_log2, S_log2), sp.Rational(-211,144)),
    sym_scale(sh(S_log2, S_log2, S_log2, S_log3), sp.Rational(101,24)),
    sym_scale(sh(S_log2, S_log2, S_log3, S_log3), sp.Rational(-15,8)),
    sym_scale(sh(S_log2, S_log3, S_log3, S_log3), sp.Rational(-35,24)),
    sym_scale(sh(S_log3, S_log3, S_log3, S_log3), sp.Rational(113,192)))
S_P = sym_expanded(S_P)
print("=== S_five + S_P:")
tot = sym_add(S_five, S_P)
for k, v in sorted(tot.items(), key=str): print("   ", k, v)
print("=== S_om:")
for k, v in sorted(S_om.items(), key=str): print("   ", k, v)
print("=== diff T:")
T = sym_sub(S_om, tot)
for k, v in sorted(T.items(), key=str): print("   ", k, v)
print()
print("S(Li_3(1/3)) =", S_Li3_13)
print("S(Li_3(-1/2)) =", S_Li3_m12)
print("S(Li_2(1/3)) =", S_Li2_13)
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
