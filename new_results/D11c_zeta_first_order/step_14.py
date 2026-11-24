import mpmath as mp
mp.mp.dps = 70

z3 = mp.zeta(3)
pi = mp.pi

# (i) P'''(1) via analytic third derivative of the arclength integral
# P(p) = 4(1/p) ∫_0^1 sqrt((1-t)^{2u} + t^{2u}) dt, u = 1/p - 1
# P'''(1)/sqrt(2) = ∫_0^1 [ -1/2 ln^3 S - 3 lnS l^2 - 9 ln^2 S - 18 l^2 - 24 lnS ] dt, S=t(1-t), l=ln(t/(1-t))
def integrand3(t):
    S = t*(1-t); l = mp.log(t/(1-t)); lnS = mp.log(S)
    return (-mp.mpf('0.5')*lnS**3 - 3*lnS*l**2 - 9*lnS**2 - 18*l**2 - 24*lnS)
I3 = mp.quad(integrand3, [0, mp.mpf('0.5'), 1])
Pppp_over_sqrt2 = I3
z3_from_perim = pi**2/3 + Pppp_over_sqrt2/6
print("zeta(3) via P'''(1)/(6 sqrt2) + pi^2/3:")
print(mp.nstr(z3_from_perim, 55))
print("error:", mp.nstr(abs(z3_from_perim - z3), 8))

# (ii) simplified J form: zeta(3) = 3 - pi^2/6 + J/2
J = mp.quad(lambda t: mp.log(t)**2 * mp.log(1-t), [0, 1])
z3_J = 3 - pi**2/6 + J/2
print("\nzeta(3) via 3 - pi^2/6 + (1/2)∫ln^2 t ln(1-t) dt:")
print(mp.nstr(z3_J, 55))
print("error:", mp.nstr(abs(z3_J - z3), 8))

# (iii) J as series: J = -2 sum_{k>=1} 1/(k(k+1)^3)
Jser = mp.nsum(lambda k: 1/(k*(k+1)**3), [1, mp.inf])
print("\nJ integral =", mp.nstr(J, 30))
print("-2 sum 1/(k(k+1)^3) =", mp.nstr(-2*Jser, 30))

# (iv) classical reduction: J = 2 zeta(3) - 6 + pi^2/3
print("\n2 zeta(3) - 6 + pi^2/3 =", mp.nstr(2*z3 - 6 + pi**2/3, 30))

# (v) independent direct check: P(p) by quadrature, 9-point FD for P'''(1)
def perim(p):
    def f(th):
        c, s = mp.cos(th), mp.sin(th)
        dxd = (2/p)*c**(2/p-1)*(-s); dyd = (2/p)*s**(2/p-1)*c
        return mp.sqrt(dxd**2+dyd**2)
    return 4*mp.quad(f, [0, mp.pi/2])
h = mp.mpf('0.003')
# 9-point stencil for f'''(1): coefficients for f(1+k h), k=-4..4
import sympy as sp
x, xp = sp.symbols('x xp')
coefs = {}
N = 4
# solve: sum_k c_k (kh)^m = delta_{m,3} 3!
import numpy as np
A = np.array([[ (k)**m for k in range(-N,N+1)] for m in range(2*N+1)], dtype=object)
bb = [0]* (2*N+1); bb[3] = 6
As = sp.Matrix(2*N+1, 2*N+1, lambda i,j: sp.Integer((i-N)**j if not(i==0 and j==0) else (1 if j==0 else 0)))
bs = sp.Matrix(bb)
cs = As.inv()*bs
fd = sum(sp.Rational(cs[i+N]) * perim(1 + i*h) for i in range(-N, N+1) if cs[i+N]!=0)/h**3
print("\nP'''(1) by 9-point FD of direct arclength quadrature:", mp.nstr(fd, 25))
print("predicted:", mp.nstr(6*mp.sqrt(2)*(z3-pi**2/3), 25))
