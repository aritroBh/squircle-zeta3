# Planner (W4): elementary proof of F23 - F32 = G0 (all u > 0) by Stokes on the triangle T = {s,t >= 0, s+t <= 1}:
#   u2u3 A - u1u4 B = d/ds[ s^u2 (1-s)^u1 * u3 t^(u3-1) (1-t)^u4 ] + d/dt[ t^u3 (1-t)^u4 * u1 s^u2 (1-s)^(u1-1) ],
#   A = s^(u2-1) t^(u3-1) (1-s)^u1 (1-t)^u4,  B = s^u2 (1-s)^(u1-1) t^u3 (1-t)^(u4-1);
#   boundary s+t = 1 gives u3 B(u1+u3, 1+u2+u4) + u1 B(1+u2+u4, u1+u3) = G(1+u1+u3)G(1+u2+u4)/G(1+u1+u2+u3+u4) = G0.
# Checks: (1) the divergence identity symbolically; (2) F23 - F32 = G0 numerically at generic u (F23 by 3F2, F32 by
# quadrature with the inner t-integral as an incomplete Beta).
import sympy as sp, mpmath as mp
s, t, u1, u2, u3, u4 = sp.symbols('s t u1 u2 u3 u4', positive=True)
A = s**(u2-1)*t**(u3-1)*(1-s)**u1*(1-t)**u4
B = s**u2*(1-s)**(u1-1)*t**u3*(1-t)**(u4-1)
div = sp.diff(s**u2*(1-s)**u1*u3*t**(u3-1)*(1-t)**u4, s) + sp.diff(t**u3*(1-t)**u4*u1*s**u2*(1-s)**(u1-1), t)
import random
random.seed(7)
expr = div - (u2*u3*A - u1*u4*B)
vals = [abs(sp.N(expr.subs({s: sp.Rational(random.randint(1, 97), 100), t: sp.Rational(random.randint(1, 97), 100),
        u1: sp.Rational(random.randint(1, 300), 100), u2: sp.Rational(random.randint(1, 300), 100),
        u3: sp.Rational(random.randint(1, 300), 100), u4: sp.Rational(random.randint(1, 300), 100)}), 50)) for _ in range(20)]
# (sympy's simplify does not close the power expressions; by hand the two cross terms -/+ u1u3 s^u2 (1-s)^(u1-1) t^(u3-1) (1-t)^u4 cancel)
print('(1) divergence identity: max |residual| at 20 random rational points (50 digits) =', max(vals))
mp.mp.dps = 30
G = mp.gamma
for u in [(0.3, 0.5, 0.7, 0.4), (1.3, 0.2, 0.9, 2.1), (0.6, 1.7, 0.25, 0.8)]:
    v1, v2, v3, v4 = [mp.mpf(x) for x in u]
    F23 = G(1+v2)*G(1+v1+v3)/G(1+v1+v2+v3)*mp.hyp3f2(-v4, v3, 1+v1+v3, 1+v3, 1+v1+v2+v3, 1)
    F32 = v1*v4*mp.quad(lambda x: x**v2*(1-x)**(v1-1)*mp.betainc(v3+1, v4, 0, 1-x), [0, 0.5, 1])
    G0 = G(1+v1+v3)*G(1+v2+v4)/G(1+v1+v2+v3+v4)
    print(f'(2) u={u}: F23 - F32 - G0 = {mp.nstr(F23 - F32 - G0, 3)}   (G0 = {mp.nstr(G0, 12)})')
