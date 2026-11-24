t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp, json, re, itertools
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def sym_Omega(c0):
    S = {}
    for w, cf in TERMS:
        S = sym_add(S, sym_scale(sym_H(w, c0), cf))
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
    return P.T.rank() == P.T.row_join(dv).rank(), P.rank()
def sym_Li4(z): return sym_scale(sym_G((0,0,0,1), z), -1)

c0 = sp.Rational(1, 3)
S_om = sym_expanded(sym_Omega(c0))
S_rhs = sym_expanded(sym_add(
    sym_scale(sym_Li4(-8), sp.Rational(-11,81)),
    sym_scale(sym_Li4(-3), sp.Rational(-2)),
    sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(20,3))))
diff = sym_sub(S_om, S_rhs)
print("S(Om(1/3)) - S(combo): terms", len(diff))
ok, rk = decomposable(diff, [sp.Integer(2), sp.Integer(3)])
print("decomposable:", ok)
print("S(Om(1/3)):")
for k, v in sorted(S_om.items(), key=str): print("   ", k, v)
print("S(combo):")
for k, v in sorted(S_rhs.items(), key=str): print("   ", k, v)
'''
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
