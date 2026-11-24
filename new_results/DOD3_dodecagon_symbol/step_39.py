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

# ---- Xi_p, Xi_q ----
Xid = pickle.load(open(os.path.join(d2,'Xi_verified_620.pkl'),'rb'))
Xip = mp.mpf(Xid['XiD_p']); Xiq = mp.mpf(Xid['XiD_q']); Xi = mp.mpf(Xid['Xi_from_S'])
print("check Xi = Xip + sqrt3 Xiq:", mp.nstr(abs(Xi - (Xip + r3*Xiq)),3))

# ---- assemble total Li3 content: 8 from Xi36 + 12 from R, merge ----
sol = pickle.load(open('xi36_li3_basis_solution.pkl','rb'))
li3tot = {}
for z,c0,c1 in zip(sol['zs'], sol['sol0'], sol['sol1']):
    li3tot[z] = (Fr(c0), Fr(c1))
R3c = pickle.load(open('R_li3_content.pkl','rb'))['li3_content']
for pt,c in R3c.items():
    p_,q_ = psplit(c)
    if pt in li3tot:
        li3tot[pt] = (li3tot[pt][0]+p_, li3tot[pt][1]+q_)
    else:
        li3tot[pt] = (p_,q_)
print("total distinct Li3 arguments:", len(li3tot))
for z in sorted(li3tot, key=pfloat):
    a0,a1 = li3tot[z]
    print("   Li3(%-14s ~%+.8f)  %s + %s sqrt3"%(str(z), pfloat(z), a0, a1))

# ---- numeric Li3 values (Re polylog) ----
Li3R = {z: mp.re(mp.polylog(3, mp.mpc(mpf_of(z)))) for z in li3tot}
V0 = Xip - sum(mp.mpf(a0.numerator)/a0.denominator*Li3R[z] for z,(a0,a1) in li3tot.items())
V1 = Xiq - sum(mp.mpf(a1.numerator)/a1.denominator*Li3R[z] for z,(a0,a1) in li3tot.items())
print("\nV0 =", mp.nstr(V0,25))
print("V1 =", mp.nstr(V1,25))

# ---- product basis (23) ----
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

for tag,V in [('V0',V0),('V1',V1)]:
    rel = mp.pslq([V]+basis, tol=mp.mpf('1e-570'), maxcoeff=10**12, maxsteps=200000)
    if rel is None:
        print(tag, "PSLQ: NO relation")
    else:
        resid = abs(rel[0]*V + sum(rel[i+1]*basis[i] for i in range(len(basis))))
        print(tag, "relation: resid =", mp.nstr(resid,3), " maxc =", max(abs(r) for r in rel), " V-coef:", rel[0])
        nz = [(names[i], rel[i+1]) for i in range(len(basis)) if rel[i+1]]
        print("   ", nz)
pickle.dump({'li3tot':li3tot,'V0':V0,'V1':V1}, open('xi_full_li3.pkl','wb'))
