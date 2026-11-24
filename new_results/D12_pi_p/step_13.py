
import mpmath as mp
import pickle
from fractions import Fraction
mp.mp.dps = 90
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi = mp.pi; L2 = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
Li4 = mp.polylog(4, mp.mpf('0.5')); Li5 = mp.polylog(5, mp.mpf('0.5'))

basis4 = [mp.mpf(1), L2, L2**2, L2**3, L2**4, pi**2, pi**2*L2, pi**2*L2**2, pi**4, z3, z3*L2, Li4]
names4 = ['1','log2','log2^2','log2^3','log2^4','pi^2','pi^2 log2','pi^2 log2^2','pi^4','zeta3','zeta3 log2','Li4(1/2)']
v = [N[4]]+basis4
rel = mp.pslq(v, tol=mp.mpf('1e-80'), maxcoeff=10**8, maxsteps=10000)
print("integer relation (coeff of N4 first):", rel)
from math import gcd
g = 0
for x in rel: g = gcd(g, x)
rel = [x//g for x in rel]
print("reduced:", rel)
print("\nN_4 =")
for i in range(len(basis4)):
    if rel[i+1]:
        c = Fraction(-rel[i+1], rel[0])
        print(f"   {c} * {names4[i]}")
val = -sum(rel[i+1]*basis4[i] for i in range(len(basis4)))/rel[0]
print("residual at 90 digits:", mp.nstr(abs(val-N[4]),3))
