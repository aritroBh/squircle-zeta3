# Integrate sc_n = s_n/sqrt2 over [0,1/2] (90 dps), compose u(d) = -d/(1+d), compare a_n to exact forms
import sympy as sp, json
from mpmath import mp, mpf, quad, log, sqrt
u, lt, l1 = sp.symbols('u lt l1')
NN = 10
def T(poly, N):
    p = sp.expand(poly); return sum(p.coeff(u, n)*u**n for n in range(0, N+1))
ex = sp.exp(u*lt).series(u, 0, NN+1).removeO().expand()
ey = sp.exp(u*l1).series(u, 0, NN+1).removeO().expand()
X = T((1+u)*ex, NN); Y = T((1+u)*ey, NN)
W = T(X*X + Y*Y, NN)
Sc = [sp.sqrt(W.coeff(u, 0))]
for n in range(1, NN+1):
    Sc.append(sp.expand((W.coeff(u, n) - sum(Sc[i]*Sc[n-i] for i in range(1, n)))/(2*Sc[0])))
scr = [sp.expand(sc/sp.sqrt(2)) for sc in Sc]     # strip overall sqrt(2): rational polys in lt, l1
mp.dps = 90
I = []
for n in range(0, NN+1):
    f = sp.lambdify((lt, l1), scr[n], 'mpmath')
    g = lambda t: f(log(t), log(1-t))
    val = quad(g, [0, mpf(1)/32, mpf(1)/4, mpf(1)/2])
    I.append(val)
# a_n = 8 * sum_{m<=n} I_m * (-1)^m C(n-1, m-1)   (for n>=1); a_0 = 8 I_0
from math import comb
a_dir = [8*I[0]]
for n in range(1, NN+1):
    a_dir.append(8*sum(I[m]*((-1)**m)*comb(n-1, m-1) for m in range(1, n+1)))
num = json.load(open('step2_num.json'))
print("INDEPENDENT arclength-integrand-expansion check (no Beta/Gamma operator formula):")
print("n : a_n (direct integrand expansion)        rel diff vs exact closed forms")
for n in range(0, NN+1):
    an = mpf(num['a'][n]); rd = abs(a_dir[n]-an)/max(1, abs(an))
    print(f"{n:2d}  {mp.nstr(a_dir[n], 45)}   {mp.nstr(rd, 3)}")
