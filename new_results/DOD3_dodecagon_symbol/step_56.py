import sys, os, pickle
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
L3c = (mp.zeta(3, mp.mpf(1)/12) - mp.zeta(3, mp.mpf(5)/12) - mp.zeta(3, mp.mpf(7)/12) + mp.zeta(3, mp.mpf(11)/12)) / 12**3
Q = Quotient()
dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
pd = pickle.load(open('prod_sym_data.pkl','rb'))
prod_syms = pd['prod_syms']
PM = sp.Matrix.hstack(*[Q.coord(s_).T for s_ in prod_syms])
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
prod_meta = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b; prod_meta.append(('log',(a,b,c)))
for w in b2zs:
    for p in ('2','3','e'): prod_meta.append(('li2',w,p))
def prodval(mta):
    if mta[0]=='log':
        a,b,c = mta[1]; return lp2**a*lp3**b*le**c
    _,w,p = mta
    return mp.polylog(2, mpf_of(w))*{'2':lp2,'3':lp3,'e':le}[p]
DM_full = sp.Matrix.hstack(*[d.T for d in Q.D])
M = sp.zeros(27,0)
for j in range(DM_full.cols):
    M2 = sp.Matrix.hstack(M, DM_full[:,j]) if M.cols else DM_full[:,j]
    if M.cols==0 or M2.rank()>M.rank(): M=M2
    if M.rank()==19: break
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

# test: T = Li3(z0) for several non-basis z0
testzs = [(3,-1,4),(1,0,4),(9,-3,4),(2,0,3),(3,2,9),(-5,3,1)]
for z0 in testzs:
    z0 = P(*z0)
    T_sym = S_Li3(z0)
    solT = BM.gauss_jordan_solve(pi(Q.coord(T_sym)))[0]
    diff = Q.coord(T_sym).T
    for j,i in enumerate(bidx):
        diff -= sp.Rational(solT[j])*gens[i].T
    solP = PM.gauss_jordan_solve(diff)[0]
    assert PM*solP == diff
    val = mp.re(mp.polylog(3, mp.mpc(mpf_of(z0))))
    for j,i in enumerate(bidx):
        if solT[j]:
            c = sp.Rational(solT[j])
            val -= mp.mpf(c.numerator)/c.denominator * mp.re(mp.polylog(3, mp.mpc(mpf_of(zs[i]))))
    for j in range(19):
        if solP[j]:
            c = sp.Rational(solP[j])
            val -= mp.mpf(c.numerator)/c.denominator * prodval(prod_meta[j])
    rel = mp.pslq([val, L3c, z3c, z2*l2, z2*l3, z2*lu], tol=mp.mpf('1e-585'), maxcoeff=10**10, maxsteps=10000)
    resid = None
    if rel:
        resid = abs(rel[0]*val + rel[1]*L3c + rel[2]*z3c + rel[3]*z2*l2 + rel[4]*z2*l3 + rel[5]*z2*lu)
    print("Li3(%s): rem=%s pslq=%s resid=%s"%(str(z0), mp.nstr(val,10), rel, mp.nstr(resid,3) if resid else None))
