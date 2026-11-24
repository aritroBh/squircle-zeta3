import sys, os, pickle
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *
Q = Quotient()
dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
pd = pickle.load(open('prod_sym_data.pkl','rb'))
prod_syms = pd['prod_syms']
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
prod_meta = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b; prod_meta.append(('log',(a,b,c)))
for w in b2zs:
    for p in ('2','3','e'): prod_meta.append(('li2',w,p))
PM = sp.Matrix.hstack(*[Q.coord(s_).T for s_ in prod_syms])
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
z0 = P(1,0,4)
Tc = Q.coord(S_Li3(z0))
pT = pi(Tc)
solT = BM.gauss_jordan_solve(pT)[0]
print("BM*solT == pT:", BM*solT == pT)
diff = Tc.T
for j,i in enumerate(bidx):
    diff -= sp.Rational(solT[j])*gens[i].T
print("pi(diff) == 0:", pi(diff.T) == sp.zeros(8,1))
solP = PM.gauss_jordan_solve(diff)[0]
print("PM*solP == diff:", PM*solP == diff)
# full tensor residual again:
resid = Tc.T
for j,i in enumerate(bidx):
    resid -= sp.Rational(solT[j])*gens[i].T
resid -= PM*solP
print("full tensor residual zero:", resid == sp.zeros(27,1))
# Hmm: print the identity component-wise for a few tensor slots to see WHERE the value mismatch could hide
# Check the claimed identity against the weight-2-derived one: Li3(1/4)... 
# Let me print solT again and the gens used, and MANUALLY recompute pi(T) coordinates
print("solT:", [str(sp.Rational(x)) for x in solT])
print("pi(T):", [str(sp.Rational(x)) for x in pT])
print("keep positions:", keep, "-> letters:", [Q.basis3[k] for k in keep])
