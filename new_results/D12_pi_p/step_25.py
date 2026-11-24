
import mpmath as mp
import pickle
from fractions import Fraction
from math import gcd
mp.mp.dps = 100
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3); z5=mp.zeta(5)
Li4=mp.polylog(4,mp.mpf('0.5')); Li5=mp.polylog(5,mp.mpf('0.5'))
PSI = pickle.load(open('/tmp/psi.pkl','rb'))['PSI']

print("exact decomposition of N_5:")
for m in range(0,5):
    for (i,j,k),v in sorted(PSI[m].items()):
        print(f"   {v} * (log2)^{k} J_{{{j},{i}}}")

W5 = [L2**5, pi**2*L2**3, pi**4*L2, z3*L2**2, z3*pi**2, z5, Li4*L2, Li5]
n5 = ['log2^5','pi^2 log2^3','pi^4 log2','z3 log2^2','z3 pi^2','z5','Li4 log2','Li5(1/2)']
def J(m,k): return mp.quad(lambda s: s**k*mp.log(1+mp.e**(-s))**m,[0,mp.inf])
for (j,i) in [(4,0),(3,1),(2,2),(1,3)]:
    t = J(j,i)
    if j==1:
        print(f"J_{{1,{i}}} = {i}! eta({i+2}) =", mp.nstr(t,15)); continue
    rel = mp.pslq([t]+W5, tol=mp.mpf('1e-90'), maxcoeff=10**10, maxsteps=20000)
    if rel and rel[0]!=0:
        g=0
        for x in rel: g=gcd(g,x)
        rel=[x//g for x in rel]
        val = -sum(rel[q+1]*W5[q] for q in range(len(W5)))/rel[0]
        terms=[(Fraction(-rel[q+1],rel[0]),n5[q]) for q in range(len(W5)) if rel[q+1]]
        print(f"J_{{{j},{i}}} (resid {mp.nstr(abs(val-t),2)}): {terms}")
    else:
        print(f"J_{{{j},{i}}}: NO weight-5 depth-1 relation -> genuinely new (depth>1) constant")
