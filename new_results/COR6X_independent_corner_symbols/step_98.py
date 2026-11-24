t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, json, re, itertools
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
W4 = [(parse_term(s)[0], parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
c0 = sp.Rational(1, 2)
# WITH ce-scaling
S_om = {}
for e, w, cf in W4:
    S_om = sym_add(S_om, sym_scale(sym_H(w, c0), cf * (c0 ** e)))
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
def decomposable(T, primes):
    basis = list(itertools.product(primes, repeat=4))
    idx = {m: i for i, m in enumerate(basis)}
    rows = []
    for u in (1, 2, 3):
        for m1 in itertools.product(primes, repeat=u):
            for m2 in itertools.product(primes, repeat=4 - u):
                r = [0]*len(basis)
                for w in shuffle_words(m1, m2): r[idx[w]] += 1
                rows.append(r)
    P = sp.Matrix(rows)
    dv = sp.Matrix([T.get(m, 0) for m in basis])
    return P.T.rank() == P.T.row_join(dv).rank()
T = sym_sub(S_om, sym_add(S_five, S_P))
print("WITH ce-scaling: S(Om) - S(five) - S(P): terms", len(T))
for k, v in sorted(T.items(), key=str): print("   ", k, v)
print("decomposable:", decomposable(T, [sp.Integer(2), sp.Integer(3)]))
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
