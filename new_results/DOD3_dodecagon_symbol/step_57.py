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
print("basis gens:", [str(zs[i]) for i in bidx])
z0 = P(1,0,4)
T_sym = S_Li3(z0)
solT = BM.gauss_jordan_solve(pi(Q.coord(T_sym)))[0]
print("solT:", [str(sp.Rational(x)) for x in solT])
diff = Q.coord(T_sym).T
for j,i in enumerate(bidx):
    diff -= sp.Rational(solT[j])*gens[i].T
solP = PM.gauss_jordan_solve(diff)[0]
print("solP:", [(prod_names[j] if False else prod_meta[j], str(sp.Rational(solP[j]))) for j in range(19) if solP[j]!=0])
# numeric piece by piece
val = mp.re(mp.polylog(3, mp.mpc(mpf_of(z0))))
print("Li3(1/4) =", mp.nstr(val, 20))
sub = mp.mpf(0)
for j,i in enumerate(bidx):
    if solT[j]:
        c = sp.Rational(solT[j])
        v = mp.re(mp.polylog(3, mp.mpc(mpf_of(zs[i]))))
        print("  a*Li3(%s): coef %s val %s contrib %s"%(str(zs[i]), c, mp.nstr(v,12), mp.nstr(mp.mpf(c.numerator)/c.denominator*v,12)))
        sub += mp.mpf(c.numerator)/c.denominator*v
print("Li3 combo:", mp.nstr(sub,15))
psub_ = mp.mpf(0)
for j in range(19):
    if solP[j]:
        c = sp.Rational(solP[j])
        v = prodval(prod_meta[j])
        print("  prod %s: coef %s val %s contrib %s"%(prod_meta[j], c, mp.nstr(v,12), mp.nstr(mp.mpf(c.numerator)/c.denominator*v,12)))
        psub_ += mp.mpf(c.numerator)/c.denominator*v
print("products:", mp.nstr(psub_,15))
print("val - combo - products =", mp.nstr(val - sub - psub_, 15))
