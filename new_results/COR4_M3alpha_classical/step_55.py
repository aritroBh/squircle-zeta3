import sympy as sp
import mpmath as mp
import numpy as np
mp.mp.dps = 50

# --- cos(3pi/8): d = 2cos(3pi/8). Relation: 2 = d^2 (2-d)(2+d)
dd = 2*sp.cos(3*sp.pi/8)
print("minpoly d (cos3pi/8):", sp.minimal_polynomial(dd, sp.Symbol('x')))
print("d^2(2-d)(2+d) - 2 ==", sp.simplify(dd**2*(2-dd)*(2+dd) - 2))
roots = [2*mp.cos(mp.pi/8), 2*mp.cos(3*mp.pi/8), -2*mp.cos(mp.pi/8), -2*mp.cos(3*mp.pi/8)]
els = [(lambda r: r), (lambda r: 2-r), (lambda r: 2+r)]
M = np.array([[mp.log(abs(f(r))) for f in els] for r in roots], dtype=float)
print("log-embedding matrix (d, 2-d, 2+d):")
print(np.round(M, 6))
print("det:", np.round(np.linalg.det(M[:3, :]), 6))

# --- cos(pi/12): relation (2+d)(2-d) = d^{-2}
dd2 = 2*sp.cos(sp.pi/12)
print("minpoly d (cos pi/12):", sp.minimal_polynomial(dd2, sp.Symbol('x')))
print("(2+d)(2-d)*d^2 ==", sp.simplify((2+dd2)*(2-dd2)*dd2**2))
embs = [mp.sqrt(2+mp.sqrt(3)), mp.sqrt(2-mp.sqrt(3)), -mp.sqrt(2+mp.sqrt(3)), -mp.sqrt(2-mp.sqrt(3))]
els2 = [(lambda r: mp.mpf(2)), (lambda r: r), (lambda r: 2+r)]
M2 = np.array([[mp.log(abs(f(r))) for f in els2] for r in embs], dtype=float)
print("log-embedding matrix (2, d, 2+d):")
print(np.round(M2, 6))
print("rank:", np.linalg.matrix_rank(M2, tol=1e-30))
