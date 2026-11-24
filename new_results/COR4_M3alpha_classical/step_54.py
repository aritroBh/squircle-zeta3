import sympy as sp
import mpmath as mp
mp.mp.dps = 50

# --- cos(3pi/8): d = 2cos(3pi/8), d^4-4d^2+2=0. Relation: 2 = d^2 (2-d)(2+d)
d = sp.AlgebraicNumber(sp.polys.polytools.minpoly(2*sp.cos(3*sp.pi/8), sp.Symbol('x')).as_expr(), sp.Symbol('x'))
print("minpoly d (cos3pi/8):", sp.minimal_polynomial(2*sp.cos(3*sp.pi/8), sp.Symbol('x')))
# verify 2 = d^2(4-d^2) i.e. d^2(2-d)(2+d) == 2
dd = 2*sp.cos(3*sp.pi/8)
print("d^2(2-d)(2+d) - 2 ==", sp.simplify(dd**2*(2-dd)*(2+dd) - 2))
# independence of {d, 2-d, 2+d}: log embeddings
roots = [2*mp.cos(mp.pi/8), 2*mp.cos(3*mp.pi/8), -2*mp.cos(mp.pi/8), -2*mp.cos(3*mp.pi/8)]
els = [(lambda r: r), (lambda r: 2-r), (lambda r: 2+r)]
import numpy as np
M = np.array([[mp.log(abs(f(r))) for f in els] for r in roots], dtype=float)
print("log-embedding matrix (rows=embeddings, cols=d,2-d,2+d):")
print(np.round(M, 6))
print("det:", np.round(np.linalg.det(M[:3, :]), 6))

# --- cos(pi/12): d = 2cos(pi/12), d^4-4d^2+1=0. Relation: (2+d)(2-d) = 4-d^2 = d^{-2}
dd2 = 2*sp.cos(sp.pi/12)
print("minpoly d (cos pi/12):", sp.minimal_polynomial(dd2, sp.Symbol('x')))
print("(2+d)(2-d)*d^2 ==", sp.simplify((2+dd2)*(2-dd2)*dd2**2))
# independence of {2, d, 2+d}: embeddings d -> +/- sqrt(2 +/- sqrt3)
embs = [mp.sqrt(2+mp.sqrt(3)), mp.sqrt(2-mp.sqrt(3)), -mp.sqrt(2+mp.sqrt(3)), -mp.sqrt(2-mp.sqrt(3))]
els2 = [(lambda r: mp.mpf(2)), (lambda r: r), (lambda r: 2+r)]
M2 = np.array([[mp.log(abs(f(r))) for f in els2] for r in embs], dtype=float)
print("log-embedding matrix (2, d, 2+d):")
print(np.round(M2, 6))
print("det(3x3):", np.round(np.linalg.det(M2[:3, :]), 6))
print("rank:", np.linalg.matrix_rank(M2, tol=1e-30))
