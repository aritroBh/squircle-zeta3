import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
from dod3_symbol import *

base = os.path.dirname(os.getcwd())
d1 = os.path.join(base,'DOD1_dodecagon'); d2 = os.path.join(base,'DOD2_dodecagon_finish')
def mpf_of(z):
    m,n,d = z
    return (mp.mpf(m) + n*mp.sqrt(3))/d
r3 = mp.sqrt(3)

# ---- Xi36_0, Xi36_1 numerically from Lam (atoms620) + SwG/SwD (gdswaps620) with exact coefficients ----
atoms620 = pickle.load(open(os.path.join(d1,'atoms620.pkl'),'rb'))
g = pickle.load(open(os.path.join(d2,'gdswaps620.pkl'),'rb'))
xi36 = pickle.load(open('xi36_symbol.pkl','rb'))
X0 = mp.mpf(0); X1 = mp.mpf(0)
for kind,key,xarg,yarg,(c0,c1) in xi36['terms']:
    if kind=='Lam':
        v = atoms620['Lam'][key]
    elif kind=='SwG':
        v = g['SwG'][key]
    else:
        v = g['SwD'][key]
    X0 += mp.mpf(c0.numerator)/c0.denominator * v
    X1 += mp.mpf(c1.numerator)/c1.denominator * v
print("Xi36_0 =", mp.nstr(X0,25))
print("Xi36_1 =", mp.nstr(X1,25))
# sanity: Xi36_0 + sqrt3 Xi36_1 should equal Xi36 = Xi - sum cT3a R_true
Xid = pickle.load(open(os.path.join(d2,'Xi_verified_620.pkl'),'rb'))
Xi = mp.mpf(Xid['Xi_from_S'])
Rq = pickle.load(open('R_all12_vals620.pkl','rb'))['Rq']
ac = pickle.load(open(os.path.join(d2,'atom_coefs_pq.pkl'),'rb'))
Rtot_true = mp.mpf(0)
for (i,j,k),(p,q) in [(kk[1],vv) for kk,vv in ac.items() if kk[0]=='T3a']:
    Rtot_true += (mp.mpf(p.numerator)/p.denominator + mp.mpf(q.numerator)/q.denominator*r3)*Rq[(i,j,k)]
Xi36_true = Xi - Rtot_true
print("check Xi36_0+sqrt3 Xi36_1 = Xi - Rtot_true:", mp.nstr(abs((X0+r3*X1)-Xi36_true),3))

# ---- TEST A: PSLQ Xi36_0 - sum8 a0 Li3 against 23-basis ----
sol = pickle.load(open('xi36_li3_basis_solution.pkl','rb'))
l2, l3, lu = mp.log(2), mp.log(3), mp.log(2+r3)
z2, z3c = mp.pi**2/6, mp.zeta(3)
basis = []; names = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b
        basis.append(l2**a*l3**b*lu**c); names.append("l2^%d l3^%d lu^%d"%(a,b,c))
basis += [z2*l2, z2*l3, z2*lu]; names += ["z2*l2","z2*l3","z2*lu"]
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
for w in b2zs:
    Lw = mp.polylog(2, mpf_of(w))
    for ln,lnm in [(l2,'l2'),(l3,'l3'),(lu,'lu')]:
        basis.append(Lw*ln); names.append("Li2(%s)*%s"%(str(w),lnm))
basis.append(z3c); names.append("zeta3")
for tag, XX, cs in [('Xi36_0',X0,sol['sol0']), ('Xi36_1',X1,sol['sol1'])]:
    L3v = [mp.re(mp.polylog(3, mp.mpc(mpf_of(z)))) for z in sol['zs']]
    T = XX - sum(mp.mpf(c.numerator)/c.denominator*L for c,L in zip(cs,L3v))
    rel = mp.pslq([T]+basis, tol=mp.mpf('1e-570'), maxcoeff=10**12, maxsteps=200000)
    if rel is None:
        print(tag, "PSLQ: NO relation   T=", mp.nstr(T,15))
    else:
        resid = abs(rel[0]*T + sum(rel[i+1]*basis[i] for i in range(len(basis))))
        print(tag, "relation: resid =", mp.nstr(resid,3), " maxc =", max(abs(r) for r in rel), " T-coef:", rel[0])
        print("   ", [(names[i], rel[i+1]) for i in range(len(basis)) if rel[i+1]])
