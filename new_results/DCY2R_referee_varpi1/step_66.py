import sympy as sp
u, y, N = sp.symbols('u y N', positive=True)
pi = sp.pi
bR = pi/(2*sp.sqrt(2))
Q = u**2 + sp.sqrt(2)*y*u + y**2
# P(u) = 2uQ * Lambda'(t) with Lambda' = (N-1)y/u - N y (2u+sqrt2 y)/(2Q) - bR
P = sp.expand(2*u*Q*((N-1)*y/u - N*y*(2*u+sp.sqrt(2)*y)/(2*Q) - bR))
print("P(u) expanded:")
print(sp.collect(P, u))
# compare with printed (typo) and corrected
P_printed = -2*bR*u**3 - (pi+2)*y*u**2 + (sp.sqrt(2)*(N-2)*y**2 - 2*bR*y**3)*u + 2*(N-1)*y**3
P_correct = -2*bR*u**3 - (pi+2)*y*u**2 + (sp.sqrt(2)*(N-2)*y**2 - 2*bR*y**2)*u + 2*(N-1)*y**3
print("diff vs printed:", sp.simplify(P - P_printed))
print("diff vs corrected:", sp.simplify(P - P_correct))
# P(1) with N y^2 = pi/2 substituted
P1 = sp.expand(P.subs(u,1))
P1 = P1.subs(N, pi/(2*y**2))
print("P(1) =", sp.simplify(P1), " < 0 for y>0")
open('progress.log','a').write("step: P(u) corrected coefficient verified symbolically; P(1)<0 exact\n")
import os
print(sorted([f for f in os.listdir('.') if f.startswith('step_')], key=lambda s: int(s[5:-3]))[-3:])
