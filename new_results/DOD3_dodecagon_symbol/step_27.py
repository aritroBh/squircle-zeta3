import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 620
from dod3_symbol import *

base = os.path.dirname(os.getcwd())
d1 = os.path.join(base,'DOD1_dodecagon'); d2 = os.path.join(base,'DOD2_dodecagon_finish')

def mpf_of(z):
    m,n,d = z
    return (mp.mpf(m) + n*mp.sqrt(3))/d

# ---- Xi (620 digits) ----
Xid = pickle.load(open(os.path.join(d2,'Xi_verified_620.pkl'),'rb'))
Xi = mp.mpf(Xid['Xi_from_S'])
print("Xi =", mp.nstr(Xi, 30))

# ---- 8 Li3 terms ----
sol = pickle.load(open('xi36_li3_basis_solution.pkl','rb'))
z8 = sol['zs']; s0 = sol['sol0']; s1 = sol['sol1']
L3v = [mp.polylog(3, mpf_of(z)) for z in z8]
r3 = mp.sqrt(3)
comb = sum((mp.mpf(int(a.numerator))/int(a.denominator) + mp.mpf(int(b.numerator))/int(b.denominator)*r3)*L
           for a,b,L in zip(s0,s1,L3v))
print("Li3-combo =", mp.nstr(comb, 30))
D = Xi - comb
print("D = Xi - combo =", mp.nstr(D, 30))

# ---- product basis ----
l2, l3, lu = mp.log(2), mp.log(3), mp.log(2+r3)
z2, z3c = mp.pi**2/6, mp.zeta(3)
import itertools
logcubes = []
names = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b
        logcubes.append(l2**a * l3**b * lu**c)
        names.append("l2^%d l3^%d lu^%d"%(a,b,c))
prods = list(logcubes)
prods += [z2*l2, z2*l3, z2*lu]; names += ["z2*l2","z2*l3","z2*lu"]
b2 = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
L2v = [mp.polylog(2, mpf_of(z)) for z in b2]
for w,L in zip(b2,L2v):
    for ln,lnm in [(l2,'l2'),(l3,'l3'),(lu,'lu')]:
        prods.append(L*ln); names.append("Li2(%s)*%s"%(str(w),lnm))
basis = prods + [z3c]; names += ["zeta3"]
print("product+const basis size:", len(basis))

# ---- PSLQ ----
rel = mp.pslq([D]+basis, tol=mp.mpf('1e-580'), maxcoeff=10**12, maxsteps=100000)
if rel is None:
    print("PSLQ: NO relation found")
else:
    resid = abs(rel[0]*D + sum(rel[i+1]*basis[i] for i in range(len(basis))))
    print("PSLQ relation found, resid =", mp.nstr(resid,3), " maxcoeff =", max(abs(r) for r in rel))
    print("D coefficient:", rel[0])
    for i,r in enumerate(rel[1:]):
        if r: print("   %+8s  %s"%(r, names[i]))
pickle.dump({'D':D,'comb':comb}, open('xi_remainder1.pkl','wb'))
