import sympy as sp
from mpmath import mp, mpf, log, sqrt, exp
mp.dps = 50
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
# numeric test at t = 0.25, u = 0.001
tv = mpf('0.25'); uv = mpf('0.001')
Sser = sum(mpf(str(Sc[n].subs({lt: log(tv), l1: log(1-tv)}).evalf(40)))*uv**n for n in range(NN+1))
ev = 1 + uv
Strue = sqrt((ev*tv**(ev-1))**2 + (ev*(1-tv)**(ev-1))**2)
print("series:", mp.nstr(Sser, 40))
print("true  :", mp.nstr(Strue, 40))
print("diff:", mp.nstr(abs(Sser-Strue), 3))
# test s_3 alone via finite difference of 3rd derivative
h = mpf('1e-6')
fv = lambda uu: sqrt(((1+uu)*tv**uu)**2 + ((1+uu)*(1-tv)**uu)**2)
d3 = (fv(2*h)-2*fv(h)+2*fv(-h)-fv(-2*h))/(2*h**3)
print("s_3 numeric:", mp.nstr(d3/6, 30), " s_3 symbolic:", mp.nstr(mpf(str(Sc[3].subs({lt: log(tv), l1: log(1-tv)}).evalf(40))), 30))
