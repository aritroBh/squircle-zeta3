import sys, os, pickle
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
from dod3_symbol import *

base = os.path.dirname(os.getcwd())
def mpf_of(z):
    m,n,d = z
    return (mp.mpf(m) + n*mp.sqrt(3))/d
r3 = mp.sqrt(3)
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}

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

# --- self-test T = Li21(x5, x6/x5) ---
T_sym = S_Li21(Xs[5], pdiv(Xs[6],Xs[5]))
solT = BM.gauss_jordan_solve(pi(Q.coord(T_sym)))[0]
diff = Q.coord(T_sym).T
for j,i in enumerate(bidx):
    diff -= sp.Rational(solT[j])*gens[i].T
solP = PM.gauss_jordan_solve(diff)[0]
assert PM*solP == diff
print("self-test exact decomposition found.")
print("Li3 coeffs:", [(str(zs[i]), str(sp.Rational(solT[j]))) for j,i in enumerate(bidx) if solT[j]!=0])
print("product coeffs:", [(prod_names[j], str(sp.Rational(solP[j]))) for j in range(19) if solP[j]!=0])

# numeric
def Li21ser(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,300000):
        xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>50 and abs(xp*ps)<mp.mpf(10)**(-630): break
    return tot
l2,l3,lu = mp.log(2),mp.log(3),mp.log(2+r3); z2,z3c = mp.pi**2/6, mp.zeta(3)
def prodval(nm):
    if nm.startswith('l'):
        a,b,c = (int(t) for t in nm.replace('l2^','').replace('l3^','').replace('lu^','').split())
        return l2**a*l3**b*lu**c
    if nm.startswith('Li2'):
        wstr, pl = nm[4:-3], nm[-1]
        w = eval(wstr)
        Lw = mp.polylog(2, mpf_of(w))
        return Lw*{'2':l2,'3':l3,'e':lu}[pl]
Tnum = Li21ser(mpf_of(Xs[5]), mpf_of(pdiv(Xs[6],Xs[5])))
val = Tnum
for j,i in enumerate(bidx):
    if solT[j]:
        c = sp.Rational(solT[j])
        val -= mp.mpf(c.numerator)/c.denominator * mp.re(mp.polylog(3, mp.mpc(mpf_of(zs[i]))))
for j in range(19):
    if solP[j]:
        c = sp.Rational(solP[j])
        val -= mp.mpf(c.numerator)/c.denominator * prodval(prod_names[j])
print("\nT - Li3combo - products =", mp.nstr(val,20))
# identify in <zeta3, z2*l2, z2*l3, z2*lu>
rel = mp.pslq([val, z3c, z2*l2, z2*l3, z2*lu], tol=mp.mpf('1e-590'), maxcoeff=10**14, maxsteps=10000)
print("invisible PSLQ:", rel)
if rel:
    resid = abs(rel[0]*val + rel[1]*z3c + rel[2]*z2*l2 + rel[3]*z2*l3 + rel[4]*z2*lu)
    print("resid:", mp.nstr(resid,3))
