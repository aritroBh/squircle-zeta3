import sys, os, pickle, itertools
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
from dod3_symbol import *
def mpf_of(z):
    m,n,d = z
    return (mp.mpf(m) + n*mp.sqrt(3))/d
r3 = mp.sqrt(3)
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+r3); z2,z3c = mp.pi**2/6, mp.zeta(3)
lp2, lp3, le = (l2+lu)/2, l3/2, lu
L2c = (mp.zeta(2, mp.mpf(1)/12) - mp.zeta(2, mp.mpf(5)/12) - mp.zeta(2, mp.mpf(7)/12) + mp.zeta(2, mp.mpf(11)/12)) / 12**2
B3 = list(itertools.product(('2','3','e'), repeat=2))
bidx9 = {t:i for i,t in enumerate(B3)}
def coord2(S):
    v = sp.zeros(1,9)
    for t,c in S.items():
        v[bidx9[t]] += sp.Rational(c.numerator, c.denominator)
    return v
def S_Li2(z):
    zv = vecfun(z); omzv = vecfun(psub(ONE,z))
    out = {}
    for a,ca in omzv.items():
        for b,cb in zv.items():
            out[(a,b)] = out.get((a,b),Fr(0)) - ca*cb
    return {k:x for k,x in out.items() if x!=0}
D2 = []
pqn = []
for p in ('2','3','e'):
    for q in ('2','3','e'):
        S = {(p,q):Fr(1),(q,p):Fr(1)}
        if p==q: S = {(p,q):Fr(2)}
        D2.append(coord2(S)); pqn.append((p,q))
DM2f = sp.Matrix.hstack(*[d.T for d in D2])
M = sp.zeros(9,0); indep=[]
for j in range(9):
    M2 = sp.Matrix.hstack(M, DM2f[:,j]) if M.cols else DM2f[:,j]
    if M.cols==0 or M2.rank()>M.rank(): indep.append(j); M=M2
    if M.rank()==6: break
DM2 = M
keep=[]
for i in range(9):
    e = sp.zeros(9,1); e[i]=1
    M2 = sp.Matrix.hstack(M,e)
    if M2.rank()>M.rank(): keep.append(i); M=M2
    if len(keep)==3: break
Ainv = M.inv()
def pi2(v): return sp.Matrix((Ainv*v.T)[6:,:])
sols = pickle.load(open('solsF_B8.pkl','rb'))
zs = sorted(sols.keys(), key=pfloat)
g2 = [coord2(S_Li2(z)) for z in zs]
Gp2 = [pi2(g) for g in g2]
def niceness(i):
    z = zs[i]; _,v2,v3,k = factor_sunit(z); fz=pfloat(z)
    return (0 if 0<fz<1 else (1 if fz<0 else 2), abs(v2)+abs(v3)+abs(k))
order = sorted(range(len(zs)), key=niceness)
bi=[]; BM=sp.zeros(3,0)
for i in order:
    M2 = sp.Matrix.hstack(BM,Gp2[i]) if BM.cols else Gp2[i]
    if BM.cols==0 or M2.rank()>BM.rank(): bi.append(i); BM=M2
    if len(bi)==3: break
z0 = P(1,0,4)
solT = BM.gauss_jordan_solve(pi2(coord2(S_Li2(z0))))[0]
diff = coord2(S_Li2(z0)).T
for j,i in enumerate(bi):
    diff -= sp.Rational(solT[j])*g2[i].T
solP = DM2.gauss_jordan_solve(diff)[0]
print("tensor check:", DM2*solP == diff)
lp = {'2':lp2,'3':lp3,'e':le}
val = mp.polylog(2, mpf_of(z0))
print("Li2 basis:", [str(zs[i]) for i in bi])
for j,i in enumerate(bi):
    if solT[j]:
        c = sp.Rational(solT[j]); print("   Li2(%s) coef %s"%(str(zs[i]), c))
        val -= mp.mpf(c.numerator)/c.denominator*mp.polylog(2, mpf_of(zs[i]))
for jj,j in enumerate(indep):
    c = sp.Rational(solP[jj])
    if c:
        p,q = pqn[j]
        print("   lp%s*lp%s coef %s"%(p,q,c))
        val -= mp.mpf(c.numerator)/c.denominator*lp[p]*lp[q]
print("weight-2 remainder:", mp.nstr(val, 20))
rel = mp.pslq([val, z2, L2c], tol=mp.mpf('1e-590'), maxcoeff=10**12, maxsteps=10000)
print("vs [zeta2, L2chi]:", rel)
