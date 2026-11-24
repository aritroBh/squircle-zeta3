
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, re, time
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe

C, cands = pipe.phase_enum(8, B=8)
pipe.print_cands(8, cands, C)
K = pipe.phase_kernel(8, C, cands)

# ---- M4i cross-checks ----
J = json.load(open("../M4i_w5_noL/L_relation.json"))
terms = J["terms"]
mp.mp.dps = 200
s2 = mp.sqrt(2)
m4_args = []
q_m4 = [Fraction(0)]*len(cands)
for key, val in terms.items():
    if key.startswith("Li5("):
        inner = key[4:-1]
        v = eval(inner, {"sqrt2": s2})
        # match to my candidates
        best = min(range(len(cands)), key=lambda i: abs(cands[i]['x'].emb(150)-v))
        assert abs(cands[best]['x'].emb(150)-v) < mp.mpf('1e-120'), (key, best)
        q_m4[best] = Fraction(val[1])
        m4_args.append((key, best))
print("M4i 13 args matched to my indices:", [b for _,b in m4_args])
missing = [k for k,_ in m4_args]
assert len(set(b for _,b in m4_args)) == 13

# (i) M4i Li5 vector in kernel? symbol mod Sym^5 == 0
d = 2
A = sp.Matrix([W.symbol_li5_dense(c['vx'], c['v1m'], d) for c in cands])
P,_ = W.sym5_basis(d)
Wann = W.int_left_kernel([list(col) for col in zip(*P)])
Wm = sp.Matrix(Wann)
Apr = A*Wm.T
qv = sp.Matrix([[sp.nsimplify(q) for q in q_m4]])
proj = qv*Apr
print("M4i Li5 vector: symbol in Sym^5 (proj==0):", proj == sp.zeros(1, proj.shape[1]))
# in kernel span?
Km = sp.Matrix(K)
aug = sp.Matrix(K + [[sp.nsimplify(q) for q in q_m4]])
print("M4i vector in row span of K:", aug.rank() == Km.rank())

# (ii) evaluate M4i RHS at 830 digits, compare with T = L/sqrt(2)
P = 830
mp.mp.dps = P+30
lam = mp.log(2); L = mp.log(1+mp.sqrt(2)); pi2 = mp.pi**2
def atom_val(key):
    if key == "z5": return mp.zeta(5)
    if key.startswith("Li5("):
        inner = key[4:-1]
        x = eval(inner, {"sqrt2": mp.sqrt(2)})
        return W.li5_series(x, P)
    # product of lam, L, pi2
    v = mp.mpf(1)
    for f in key.split("*"):
        v *= {"lam": lam, "L": L, "pi2": pi2}[f]
    return v
t0=time.time()
rhs = sum(Fraction(val[1])*atom_val(key) for key, val in terms.items())
Tv = W.Lchi_hurwitz(8, P)/mp.sqrt(2)
mp.mp.dps = 60
print("M4i RHS - L/sqrt2 at 830d:", mp.nstr(rhs - Tv, 5), f"({time.time()-t0:.0f}s)")
