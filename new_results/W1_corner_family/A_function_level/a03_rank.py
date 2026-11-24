"""a03_rank.py -- exact rank of the four obstruction vectors ob(F_0 o c), ob(F_0 o s), ob(F_-1 o c), ob(F_-1 o s)
in Lambda^2(Lambda^2 V_tau) (tau-alphabet t0..t3).  If rank 4, no combination with prefactors evaluated at
tau = sqrt2-1 can cancel the non-classical part, so the complementary-angle involution cannot rescue Route A."""
import sympy as sp
exec(open('a02_involution_obstruction.py').read().split("for e in sorted(F):")[0])
vecs = []
for e in (0, -1):
    S = csymbol(F[e])
    vecs.append(delta22(pull(S, CL))); vecs.append(delta22(pull(S, SL)))
keys = sorted(set().union(*vecs))
M = sp.Matrix([[sp.Rational(v.get(k, 0)) for v in vecs] for k in keys])
print("number of Lambda^2Lambda^2 coordinates touched:", len(keys))
print("rank of {ob(F0oc), ob(F0os), ob(F-1oc), ob(F-1os)} =", M.rank())
