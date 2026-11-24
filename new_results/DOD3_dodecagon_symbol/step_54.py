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
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
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

def serA(x,y):  # = sum_{m>k} y^m x^k/(m^2 k)  = Li21ser(y,x)
    tot=mp.mpf(0); ps=mp.mpf(0); yp=mp.mpf(1); xp=mp.mpf(1)
    for m in range(1,300000):
        yp*=y
        if m>=2: tot += yp*ps/m**2
        xp*=x; ps += xp/m
        if m>50 and abs(yp*ps)<mp.mpf(10)**(-630): break
    return tot
def serB(x,y):  # = sum_{m>k} y^m x^k/(m k^2)
    tot=mp.mpf(0); ps=mp.mpf(0); yp=mp.mpf(1); xp=mp.mpf(1)
    for m in range(1,300000):
        yp*=y; xp*=x
        if m>=2: tot += yp*ps/m
        ps += xp/m**2
        if m>50 and abs(yp*ps)<mp.mpf(10)**(-630): break
    return tot

xa, ya = Xs[5], pdiv(Xs[6],Xs[5])
T_sym = S_Li21(xa, ya)
solT = BM.gauss_jordan_solve(pi(Q.coord(T_sym)))[0]
diff = Q.coord(T_sym).T
for j,i in enumerate(bidx):
    diff -= sp.Rational(solT[j])*gens[i].T
solP = PM.gauss_jordan_solve(diff)[0]
assert PM*solP == diff
for tag, Tnum in [('A: S(x,y) = sum y^m x^k/(m^2 k)', serA(mpf_of(xa),mpf_of(ya))),
                  ('B: S(x,y) = sum y^m x^k/(m k^2)', serB(mpf_of(xa),mpf_of(ya)))]:
    val = Tnum
    for j,i in enumerate(bidx):
        if solT[j]:
            c = sp.Rational(solT[j])
            val -= mp.mpf(c.numerator)/c.denominator * mp.re(mp.polylog(3, mp.mpc(mpf_of(zs[i]))))
    for j in range(19):
        if solP[j]:
            c = sp.Rational(solP[j])
            val -= mp.mpf(c.numerator)/c.denominator * prodval(prod_meta[j])
    rel = mp.pslq([val, L3c, z3c, z2*l2, z2*l3, z2*lu], tol=mp.mpf('1e-585'), maxcoeff=10**10, maxsteps=10000)
    print(tag, " T=%s"%mp.nstr(Tnum,12))
    print("   rem=%s pslq=%s"%(mp.nstr(val,10), rel))
