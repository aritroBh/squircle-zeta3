# s17: exact edge data and exact edge-law coefficient for hexagon A (regular normals, weights 1,2,3).
import sympy as sp
from s04_exact_local import polygon_exact, edge_exact
r3 = sp.sqrt(3)
ns = [(0, 1), (r3/2, sp.Rational(1, 2)), (r3/2, -sp.Rational(1, 2))]   # b = 0, -pi/3, -2pi/3
ws = [1, 2, 3]
dirs = [(sp.cos(k*sp.pi/3), sp.sin(k*sp.pi/3)) for k in range(6)]
V = polygon_exact(ns, ws, [0, 0, 0], dirs)
print("vertices:", V)
Z = 0
for q in range(6):
    dat, nl = edge_exact(V[q], V[(q+1) % 6], ns, ws, [0, 0, 0])
    A, B, L, h = dat['alpha'], dat['beta'], dat['L'], dat['h']
    u = -dat['s0']/L; v = 1 - u
    Ze = sp.nsimplify(sp.simplify(6*A*B*h**2/L*(3*(A*u + B*v) - 2)))
    Z += Ze
    print(f"edge {q}: alpha={A}, beta={B}, L={sp.simplify(L)}, h={sp.simplify(h)}, u={sp.simplify(u)}, r={[sp.simplify(x['r']) for x in nl]}, Z_e={Ze}")
Z = sp.nsimplify(sp.simplify(Z))
print("Z(hexA) exact =", Z, "=", sp.N(Z, 30))
