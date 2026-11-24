import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
from dod3_symbol import *

base = os.path.dirname(os.getcwd())
d1 = os.path.join(base,'DOD1_dodecagon')
def mpf_of(z):
    m,n,d = z
    return (mp.mpf(m) + n*mp.sqrt(3))/d
r3 = mp.sqrt(3)
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}

def Li21ser(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,300000):
        xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>50 and abs(xp*ps)<mp.mpf(10)**(-630): break
    return tot

atoms620 = pickle.load(open(os.path.join(d1,'atoms620.pkl'),'rb'))
# check Lam(5,6) = Li21(x_5, x_6/x_5) convergent
v_ser = Li21ser(mpf_of(Xs[5]), mpf_of(pdiv(Xs[6],Xs[5])))
v_stored = atoms620['Lam'][(5,6)]
print("Lam(5,6): series =", mp.nstr(v_ser,20))
print("          stored =", mp.nstr(v_stored,20), " diff =", mp.nstr(abs(v_ser-v_stored),3))

# SELF-TEST: symbol identity for T = Li21(x5, x6/x5) vs Li3 span
T_sym = S_Li21(Xs[5], pdiv(Xs[6],Xs[5]))
dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
Q = Quotient()
# reuse projection from step_22 approach: rebuild quickly
DM_full = sp.Matrix.hstack(*[d.T for d in Q.D])
M = sp.zeros(27,0); indep=[]
for j in range(DM_full.cols):
    M2 = sp.Matrix.hstack(M, DM_full[:,j]) if M.cols else DM_full[:,j]
    if M.cols==0 or M2.rank()>M.rank(): indep.append(j); M=M2
    if M.rank()==19: break
D19=M
keep=[]
for i in range(27):
    e = sp.zeros(27,1); e[i]=1
    M2 = sp.Matrix.hstack(M,e)
    if M2.rank()>M.rank(): keep.append(i); M=M2
    if len(keep)==8: break
Ainv = M.inv()
def pi(v): return sp.Matrix((Ainv*v.T)[19:,:])
Gp = [pi(g) for g in gens]
def niceness(i):
    z = zs[i]; _,v2,v3,k = factor_sunit(z); fz=pfloat(z)
    return (0 if 0<fz<1 else (1 if fz<0 else 2), abs(v2)+abs(v3)+abs(k))
order = sorted(range(len(zs)), key=niceness)
bidx=[]; BM=sp.zeros(8,0)
for i in order:
    M2 = sp.Matrix.hstack(BM,Gp[i]) if BM.cols else Gp[i]
    if BM.cols==0 or M2.rank()>BM.rank(): bidx.append(i); BM=M2
    if len(bidx)==8: break
tT = pi(Q.coord(T_sym))
solT = BM.gauss_jordan_solve(tT)[0]
assert BM*solT == tT
print("self-test Li21(x5,x6/x5): coefficients:")
L3v = []
for j,i in enumerate(bidx):
    print("   z=%-14s c=%s"%(str(zs[i]), sp.Rational(solT[j])))
    L3v.append(mp.re(mp.polylog(3, mp.mpc(mpf_of(zs[i])))))
Tnum = v_ser
Trem = Tnum - sum(mp.mpf(sp.Rational(solT[j]).numerator)/sp.Rational(solT[j]).denominator*L for j,L in enumerate(L3v))
print("T - combo =", mp.nstr(Trem,15))
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+r3); z2,z3c = mp.pi**2/6, mp.zeta(3)
basis=[];names=[]
for a in range(4):
    for b in range(4-a):
        c=3-a-b; basis.append(l2**a*l3**b*lu**c); names.append("l2^%d l3^%d lu^%d"%(a,b,c))
basis += [z2*l2,z2*l3,z2*lu]; names += ["z2*l2","z2*l3","z2*lu"]
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
for w in b2zs:
    Lw = mp.polylog(2,mpf_of(w))
    for ln,lnm in [(l2,'l2'),(l3,'l3'),(lu,'lu')]:
        basis.append(Lw*ln); names.append("Li2(%s)*%s"%(str(w),lnm))
basis.append(z3c); names.append("zeta3")
rel = mp.pslq([Trem]+basis, tol=mp.mpf('1e-570'), maxcoeff=10**12, maxsteps=200000)
if rel is None: print("SELF-TEST PSLQ: NO relation")
else:
    resid = abs(rel[0]*Trem + sum(rel[i+1]*basis[i] for i in range(len(basis))))
    print("SELF-TEST relation: resid =", mp.nstr(resid,3), " maxc =", max(abs(r) for r in rel))
    print("   ", [(names[i], rel[i+1]) for i in range(len(basis)) if rel[i+1]])
