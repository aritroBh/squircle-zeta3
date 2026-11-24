# Quick sanity check of the exact decomposition pieces (independent re-derivation):
# A1 := n^3 * J where J = ∬ x^{n-1} y^{n-1} sqrt(x^2+y^2) dxdy
# claim: J = 2/(2n+1) ∫_0^1 u^{n-1} sqrt(1+u^2) du  =>  A1 = n^2 Φ(n)/(1+1/(2n))
from mpmath import mp, mpf, sqrt, nstr, quad, inf, exp
mp.dps = 40
def Phi(n):
    n = mpf(n)
    return quad(lambda y: exp(-n*y)*sqrt(1+exp(-2*y)), [0, mpf(20)/n, inf])
def J_direct(n):
    n = mpf(n)
    return quad(lambda x,y: x**(n-1)*y**(n-1)*sqrt(x*x+y*y), [0,1],[0,1])
for n in [3, 7]:
    Jd = J_direct(n)
    Jform = 2/(2*n+1)*quad(lambda u: u**(n-1)*sqrt(1+u*u), [0,1])
    print(n, 'J diff:', nstr(Jd-Jform, 8), ' n^3 J vs n^2 Φ/(1+1/2n):', nstr(n**3*Jd - n**2*Phi(n)/(1+mpf(1)/(2*n)), 8))
