
import mpmath as mp
import pickle
from fractions import Fraction
from math import gcd
mp.mp.dps = 100
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3); z5=mp.zeta(5); z7=mp.zeta(7)
Li4=mp.polylog(4,mp.mpf('0.5')); Li5=mp.polylog(5,mp.mpf('0.5')); Li6=mp.polylog(6,mp.mpf('0.5'))
Li7=mp.polylog(7,mp.mpf('0.5'))

r = {k: 4*(N[k-1]-(-L2)**k/mp.factorial(k)) for k in range(3,42)}
a = {m: sum(L2**j/mp.factorial(j)*r[m-j] for j in range(0,m-2)) for m in range(3,42)}

basis = [pi**2, z3, pi**4, z5, pi**6, z7, pi**8, z3**2, z3*z5, z3**2*pi**2, z3*pi**4,
         L2, L2**2, L2**3, L2**4, L2**5, L2**6, Li4, Li5, Li6, Li7,
         z3*L2, z5*L2, pi**2*L2, pi**4*L2, z3**2*L2, z3*pi**2*L2, pi**2*L2**2, z3*L2**2, Li4*L2, Li5*L2]
names = ['pi^2','z3','pi^4','z5','pi^6','z7','pi^8','z3^2','z3 z5','z3^2 pi^2','z3 pi^4',
         'log2','log2^2','log2^3','log2^4','log2^5','log2^6','Li4','Li5','Li6','Li7',
         'z3 log2','z5 log2','pi^2 log2','pi^4 log2','z3^2 log2','z3 pi^2 log2','pi^2 log2^2','z3 log2^2','Li4 log2','Li5 log2']

for m in [6,7,8,9,10]:
    t = a[m]-a[m-1]
    rel = mp.pslq([t]+basis, tol=mp.mpf('1e-80'), maxcoeff=10**12, maxsteps=20000)
    if rel and rel[0]!=0:
        g=0
        for x in rel: g=gcd(g,x)
        rel=[x//g for x in rel]
        terms=[(Fraction(-rel[i+1],rel[0]),names[i]) for i in range(len(basis)) if rel[i+1]]
        val = -sum(rel[i+1]*basis[i] for i in range(len(basis)))/rel[0]
        print(f"a_{m}-a_{m-1} (resid {mp.nstr(abs(val-t),2)}): ", terms)
    else:
        print(f"a_{m}-a_{m-1}: no relation")
