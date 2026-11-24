# Copied from research/08_theory/K1f_hexagon_order5/step_27.py (our notes) with the working directory taken from the environment variable G9WORK.
# Only paths and the module name of the e-sum engine were changed.
import os

# step_27.py (v2): local Beta table Jloc(N;a,b), a+b<=5, symbolic in N.
# lnGamma(m+z)-lnGamma(m) uses H_{m-1}: A-series H_N = Hs - DEN^r, C-series H_{N+1} = Hs.
import sympy as sp, pickle
WD = os.environ['G9WORK']
u, v = sp.symbols('u v')
z2, z3, z5 = sp.symbols('z2 z3 z5')
z4 = sp.Rational(2,5)*z2**2
Hs = sp.symbols('H1 H2 H3 H4 H5')          # H_{N+1}^{(r)}
T2 = sp.symbols('T1 T2 T3 T4 T5')          # 1/(N+2)^r
DEN = sp.Symbol('DEN')                     # 1/(N+1)

def lg(h):
    return [h[0], (z2 - h[1])/2, (-2*z3 + 2*h[2])/6, (6*z4 - 6*h[3])/24, (-24*z5 + 24*h[4])/120]
A = lg([Hs[r] - DEN**(r+1) for r in range(5)])     # H_N^{(r)} = H_{N+1}^{(r)} - 1/(N+1)^r
B = lg([0,0,0,0,0])                                 # H_0 = 0 (pure zeta)
C = lg(list(Hs))                                    # H_{N+1}
S = sp.expand(sum(A[j]*u**(j+1) + B[j]*v**(j+1) - C[j]*(u+v)**(j+1) for j in range(5)))
E = sp.Integer(1); P = sp.Integer(1)
for j in range(1, 6):
    P = sp.expand(P*S)
    E = sp.expand(E + P/sp.factorial(j))
tab = {}
for a in range(6):
    for b in range(6-a):
        co = sp.expand(E.coeff(u, a).coeff(v, b))
        tab[(a, b)] = sp.expand(DEN*sp.factorial(a)*sp.factorial(b)*co)
print("Jloc(N;1,0) =", tab[(1,0)], " (should be -DEN^2)")
print("Jloc(N;0,1) =", tab[(0,1)])
pickle.dump({k: str(v) for k,v in tab.items()}, open(f'{WD}/jloc_sym.pkl','wb'))
print("saved jloc_sym.pkl")

import mpmath as mp
mp.mp.dps = 50
zv = {z2: sp.pi**2/6, z3: sp.zeta(3), z5: sp.zeta(5)}
import random
random.seed(1)
worst = mp.mpf(0)
for trial in range(12):
    N = random.randint(0, 40)
    subs = {Hs[r-1]: sp.harmonic(N+1, r) for r in range(1,6)}
    subs.update({T2[r-1]: sp.Rational(1,(N+2)**r) for r in range(1,6)})
    subs[DEN] = sp.Rational(1, N+1)
    subs.update(zv)
    for (a,b), e in tab.items():
        val = mp.mpf(sp.N(e.subs(subs, simultaneous=True), 50))
        num = mp.quad(lambda t: t**N*mp.log(t)**a*mp.log(1-t)**b, [0, mp.mpf('0.5'), 1])
        d = abs(val-num); worst = max(worst, d)
print("verification: 12 random N x 21 (a,b), worst |diff| =", mp.nstr(worst, 3))
