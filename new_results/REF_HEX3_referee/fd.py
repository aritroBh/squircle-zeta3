"""Finite-difference check of P'''(1) (no differentiation under the integral): perimeters at p=1+-h,1+-2h,
5-point stencil, Richardson extrapolation in h^2."""
import sys, sympy as sp, mpmath as mp, direct, formula
from run_tests import HEX
mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 60
name = sys.argv[1]
V, w = HEX[name]
Vm = [(mp.mpf(sp.N(x, mp.mp.dps+20)), mp.mpf(sp.N(y, mp.mp.dps+20))) for x, y in V]
H = direct.setup(Vm, [mp.mpf(sp.N(x, mp.mp.dps+20)) for x in w])
ref = 2*sum(formula.E_formula(ed) for ed in formula.edges_data(V, w))
h0 = mp.mpf(1)/32; levels = 7
R = []
for k in range(levels):
    h = h0/2**k
    P = {s: direct.perimeter(H, 1+s*h) for s in (-2,-1,1,2)}
    D = (P[2]-2*P[1]+2*P[-1]-P[-2])/(2*h**3)
    row = [D]
    for j in range(1, k+1):
        row.append((4**j*row[j-1]-R[k-1][j-1])/(4**j-1))
    R.append(row)
    print(k, mp.nstr(row[-1], 40), ' minus table:', mp.nstr(row[-1]-ref, 5)); sys.stdout.flush()
