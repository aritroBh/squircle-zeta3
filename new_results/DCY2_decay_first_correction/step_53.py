import sympy as sp
u = sp.Symbol('u')
c = [1, sp.Rational(-1,2), sp.Rational(5,24), sp.Rational(-1,16), sp.Rational(79,5760)]
for k in range(1,4):
    g_k = sp.simplify(c[k]*sp.gamma(sp.Rational(k)+sp.Rational(3,2))/sp.gamma(sp.Rational(3,2)))
    print(f"g_{k} = {g_k} = {sp.N(g_k,12)}")
# C_rho for u in [0,1]: bound |v^{(4)}|/24 -- find max numerically
from mpmath import mp, mpf, sqrt, exp, diff
mp.dps = 50
def v(x):
    return sqrt((1-exp(-2*x))/(2*x)) if x != 0 else mpf(1)
# max |v''''| on [0,1] via fine grid
mx = 0
for i in range(2001):
    x = mpf(i)/2000
    d4 = diff(v, x, 4)
    mx = max(mx, abs(d4))
print("max |v4| on [0,1] ~", mx, " C_rho ~", mx/24)
# tail constant C_rho' : |rho_2(u)| for u>=1: rho_2 = v - 1 + u/2 - 5u^2/24... wait c1=-1/2,c2=5/24
# rho_2 = v - (1 - u/2 + 5u^2/24): bound |rho_2| <= C_rho2' u^2 for u >= 1?
def rho2(x):
    return v(x) - (1 - x/2 + 5*x**2/24)
vals = [(mpf(k), abs(rho2(mpf(k)))/mpf(k)**2) for k in [1,2,5,10,100]]
print(vals)
