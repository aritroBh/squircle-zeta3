
import mpmath as mp
from fractions import Fraction
mp.mp.dps = 100
pi = mp.pi; L2 = mp.log(2); z3 = mp.zeta(3)
Li4 = mp.polylog(4, mp.mpf('0.5'))

def J(m,k):
    return mp.quad(lambda s: s**k * mp.log(1+mp.e**(-s))**m, [0, mp.inf])

W4 = [L2**4, pi**2*L2**2, pi**4, z3*L2, Li4]
n4 = ['log2^4','pi^2 log2^2','pi^4','zeta3 log2','Li4(1/2)']
for (m,k) in [(3,0),(2,1)]:
    t = J(m,k)
    rel = mp.pslq([t]+W4, tol=mp.mpf('1e-90'), maxcoeff=10**10, maxsteps=10000)
    if rel:
        from math import gcd
        g=0
        for x in rel: g=gcd(g,x)
        rel=[x//g for x in rel]
        s = " + ".join(f"{Fraction(-rel[i+1],rel[0])}*{n4[i]}" for i in range(len(W4)) if rel[i+1])
        val = -sum(rel[i+1]*W4[i] for i in range(len(W4)))/rel[0]
        print(f"J_{{{m},{k}}}: residual={mp.nstr(abs(val-t),3)}")
        print(f"   J_{{{m},{k}}} = {s}")
    else:
        print(f"J_{{{m},{k}}}: no relation")
