
import mpmath as mp
import pickle
from fractions import Fraction
mp.mp.dps = 100
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3); z5=mp.zeta(5); z7=mp.zeta(7)
Li4=mp.polylog(4,mp.mpf('0.5')); Li5=mp.polylog(5,mp.mpf('0.5')); Li6=mp.polylog(6,mp.mpf('0.5'))
Li7=mp.polylog(7,mp.mpf('0.5'))

r = {k: 4*(N[k-1]-(-L2)**k/mp.factorial(k)) for k in range(3,42)}
a = {}
for m in range(3,42):
    a[m] = sum(L2**j/mp.factorial(j)*r[m-j] for j in range(0,m-2))

# basis of Euler sums weights up to 7
basis = [pi**2, z3, pi**4, z3*L2, pi**2*L2**2, L2**4, Li4, z5, pi**2*z3, z3**2, z3*pi**2, pi**4*L2,
         z5*L2, Li5, pi**6, z7, z3*z5 if False else mp.zeta(5)*L2**0, L2**5, pi**2*L2**3, z3*L2**2, Li4*L2,
         pi**4*L2**2, z3*pi**2*L2, z5*pi**2, z3**2*L2, Li5*L2, Li6, pi**2*L2**4, L2**6, z3*L2**3, pi**6*L2, z7*L2,
         Li7, pi**8, z3*pi**4, z5*pi**2*L2**0, z3**2*pi**2]
names = ['pi^2','z3','pi^4','z3 log2','pi^2 log2^2','log2^4','Li4','z5','pi^2 z3','z3^2','z3 pi^2','pi^4 log2',
         'z5 log2','Li5','pi^6','z7','z5','log2^5','pi^2 log2^3','z3 log2^2','Li4 log2',
         'pi^4 log2^2','z3 pi^2 log2','z5 pi^2','z3^2 log2','Li5 log2','Li6','pi^2 log2^4','log2^6','z3 log2^3','pi^6 log2','z7 log2',
         'Li7','pi^8','z3 pi^4','z5 pi^2','z3^2 pi^2']
# try incremental: a_m - a_{m-1}
for m in [6,7,8,9,10]:
    t = a[m]-a[m-1]
    rel = mp.pslq([t]+basis, tol=mp.mpf('1e-85'), maxcoeff=10**12, maxsteps=20000)
    if rel:
        from math import gcd
        g=0
        for x in rel: g=gcd(g,x)
        rel=[x//g for x in rel]
        terms=[(Fraction(-rel[i+1],rel[0]),names[i]) for i in range(len(basis)) if rel[i+1]]
        val = -sum(rel[i+1]*basis[i] for i in range(len(basis)))/rel[0]
        print(f"a_{m}-a_{m-1} (resid {mp.nstr(abs(val-t),2)}): ", terms)
    else:
        print(f"a_{m}-a_{m-1}: no relation found")
