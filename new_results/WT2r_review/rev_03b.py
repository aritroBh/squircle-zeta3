
import sys
import sympy as sp
u = sp.symbols('u'); I = sp.I; S2 = sp.sqrt(2)
tau10 = (-1-I)/2
a = 1+S2; b = S2-1
DL1s = 1/(u+I-b) + 1/(u+I+a) - 1/(u+I)
for dd in [2,3]:
    gs = sp.symbols('g1:'+str(2*dd))
    G = sum(gs[j-1]*u**(-j) for j in range(1,2*dd))
    prod = G * DL1s
    ser = sp.series(prod, u, 0, 1)
    print("d=",dd,"series:", ser)
    coeff = sp.expand(ser.removeO()).coeff(u, -(2*dd-2))
    print("   [u^-(2d-2)] =", sp.radsimp(sp.simplify(coeff)), " expected tau1_0*g_top =", sp.simplify(tau10*gs[2*dd-2]))
