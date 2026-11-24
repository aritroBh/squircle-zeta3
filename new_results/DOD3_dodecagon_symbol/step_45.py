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
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']

prod_meta = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b
        prod_meta.append(('log',(a,b,c)))
for w in b2zs:
    for p in ('2','3','e'):
        prod_meta.append(('li2',w,p))
pd = pickle.load(open('prod_sym_data.pkl','rb'))
prod_syms, prod_names = pd['prod_syms'], pd['prod_names']
Q = Quotient()
PM = sp.Matrix.hstack(*[Q.coord(s).T for s in prod_syms])
dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
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

T_sym = S_Li21(Xs[5], pdiv(Xs[6],Xs[5]))
solT = BM.gauss_jordan_solve(pi(Q.coord(T_sym)))[0]
diff = Q.coord(T_sym).T
for j,i in enumerate(bidx):
    diff -= sp.Rational(solT[j])*gens[i].T
solP = PM.gauss_jordan_solve(diff)[0]
assert PM*solP == diff

def Li21ser(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,300000):
        xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>50 and abs(xp*ps)<mp.mpf(10)**(-630): break
    return tot
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+r3); z2,z3c = mp.pi**2/6, mp.zeta(3)
def prodval(mta):
    if mta[0]=='log':
        a,b,c = mta[1]; return l2**a*l3**b*lu**c
    _,w,p = mta
    return mp.polylog(2, mpf_of(w))*{'2':l2,'3':l3,'e':lu}[p]
Tnum = Li21ser(mpf_of(Xs[5]), mpf_of(pdiv(Xs[6],Xs[5])))
val = Tnum
for j,i in enumerate(bidx):
    if solT[j]:
        c = sp.Rational(solT[j])
        val -= mp.mpf(c.numerator)/c.denominator * mp.re(mp.polylog(3, mp.mpc(mpf_of(zs[i]))))
for j in range(19):
    if solP[j]:
        c = sp.Rational(solP[j])
        val -= mp.mpf(c.numerator)/c.denominator * prodval(prod_meta[j])
print("T - Li3combo - products =", mp.nstr(val,20))
rel = mp.pslq([val, z3c, z2*l2, z2*l3, z2*lu], tol=mp.mpf('1e-590'), maxcoeff=10**14, maxsteps=10000)
if rel:
    resid = abs(rel[0]*val + rel[1]*z3c + rel[2]*z2*l2 + rel[3]*z2*l3 + rel[4]*z2*lu)
    print("invisible PSLQ:", rel, " resid:", mp.nstr(resid,3))
else:
    print("invisible PSLQ: none")
