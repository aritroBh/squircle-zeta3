t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, json, re, itertools
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
W4 = [(parse_term(s)[0], parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def sym_Omega_scaled(c0):
    S = {}
    for e, w, cf in W4:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf * (c0 ** e)))
    return S
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
def sym_Li4(z): return sym_scale(sym_G((0,0,0,1), z), -1)

c0 = sp.Rational(1, 3)
S_om = sym_expanded(sym_Omega_scaled(c0))
S_rhs = sym_expanded(sym_add(
    sym_scale(sym_Li4(-8), sp.Rational(-11,81)),
    sym_scale(sym_Li4(-3), sp.Rational(-2)),
    sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(20,3))))
diff = sym_sub(S_om, S_rhs)
print("1/3 with ce-scaling: S(Om) - S(combo): terms", len(diff),
      " decomposable:", decomposable(diff, [sp.Integer(2), sp.Integer(3)]))

# also 1/2 with the 3-term Lie class (planner: 3-term == 5-term mod shuffles)
c0 = sp.Rational(1, 2)
S_om2 = sym_expanded(sym_Omega_scaled(c0))
S_3t = sym_expanded(sym_add(
    sym_scale(sym_Li4(-2), sp.Rational(-21)),
    sym_scale(sym_Li4(sp.Rational(1,3)), sp.Rational(13,2)),
    sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(13,2))))
d2 = sym_sub(S_om2, S_3t)
print("1/2 with ce-scaling: S(Om) - S(-21Li4(-2)+(13/2)Li4(1/3)+(13/2)Li4(2/3)): terms", len(d2),
      " decomposable:", decomposable(d2, [sp.Integer(2), sp.Integer(3)]))
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
