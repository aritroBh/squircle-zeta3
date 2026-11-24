import sympy as sp
from fractions import Fraction as Fr
m = sp.symbols('m', integer=True, positive=True)
pi, L = sp.pi, sp.log(2)
H = sp.harmonic
def E_2m(): return H(m)-H(2*m)
def E_2m1(): return H(m-1)-H(2*m-2)+1/(2*m-1)
def Ktval(k, Ek_1):
    K1 = pi/8 - pi*L/4
    S = (pi/4)*(-1)**k*(1 - 2*Ek_1 - (-1)**(k-1)/k)
    return 2*(-1)**(k-1)*K1 + S
def B_(mm, k): return (-1)**k*pi/(2*mm) - pi*L/2
def Tv(mm, k, Ek_1): return -2*k*Ktval(k, Ek_1) - pi/2 + pi*L/2 - pi/(4*k)
v1 = (pi/2 + B_(4*m, 2*m) + Tv(4*m, 2*m, E_2m1()))/(4*m)
v2 = (pi/2 + B_(4*m+2, 2*m+1) + Tv(4*m+2, 2*m+1, E_2m()))/(4*m+2)
expr1 = (v1+v2)/pi
for mv in [1,2,3]:
    print(mv, sp.N(expr1.subs(m,mv),15), " law:", sp.N(1/(8*m*(2*m+1)**2).subs(m,mv),15))
