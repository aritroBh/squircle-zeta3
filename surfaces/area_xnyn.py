"""Area of the surface z = x^n y^n over the unit square (Ahmed's Problem 4), as n -> infinity.

Prints E(n) = A(n) - 1 and n*E(n) - sqrt(2) for n = 10 ... 10^5, which checks the leading term
    A(n) = 1 + sqrt(2)/n + O(log^2 n / n^2).
Substitution x = e^{-s}, y = e^{-t}; the inner integral runs along s + t = sigma.  mpmath, 25 digits.
Run:  python area_xnyn.py
"""
from mpmath import *
mp.dps=25
def E(n):
    n=mpf(n); e=1/n
    # u=e^{-s}, v=e^{-t}; sigma=s+t, s in [0,sigma]
    def inner(sig):
        f=lambda s: (sqrt(1+n**2*exp(-2*sig)*(exp(2*e*s)+exp(2*e*(sig-s))))-1)
        return quad(f,[0,sig/2,sig])*exp(-e*sig)
    c=log(n)
    return quad(inner,[0,c-3,c,c+3,c+10,2*c+40,inf])/n**2
for n in [10,20,100,1000,10000,100000]:
    v=E(n); print(n, v, n*v, n*v-sqrt(2))
