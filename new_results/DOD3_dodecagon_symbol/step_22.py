import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
Q = Quotient()

# ---- projection pi: Q^27 -> Q^8 complement of Decomp ----
t0=time.time()
DM = sp.Matrix.hstack(*[d.T for d in Q.D])   # 27 x 27, rank 19
keep = []
M = DM
for i in range(27):
    e = sp.zeros(27,1); e[i]=1
    M2 = sp.Matrix.hstack(M, e)
    if M2.rank() > M.rank():
        keep.append(i); M = M2
    if len(keep)==8: break
assert len(keep)==8
E = sp.zeros(27,8)
for j,i in enumerate(keep): E[i,j]=1
A = sp.Matrix.hstack(DM, E)   # 27x27 invertible
Ainv = A.inv()
print("projection built, kept coordinate positions:", keep, " t=%.1fs"%(time.time()-t0))
def pi(v):  # v: 1x27 row -> 8x1 col of complement coords
    w = Ainv * v.T
    return w[19:,:]
Gp = [pi(g) for g in gens]
# greedy basis among gens, preference: |z float| in nice range, small exponents
def niceness(i):
    z = zs[i]; _,v2,v3,k = sols_fac = factor_sunit(z)
    fz = pfloat(z)
    pref = 0 if 0 < fz < 1 else (1 if fz < 0 else 2)
    return (pref, abs(v2)+abs(v3)+abs(k))
order = sorted(range(len(zs)), key=niceness)
basis_idx = []
Bcols = []
BM = sp.zeros(8,0)
for i in order:
    M2 = sp.Matrix.hstack(BM, Gp[i]) if BM.cols else Gp[i]
    if BM.cols==0 or M2.rank()>BM.rank():
        basis_idx.append(i); BM = M2
    if len(basis_idx)==8: break
print("basis gens (8):", basis_idx, "rank:", BM.rank())
for i in basis_idx:
    print("   z = %-14s ~ %+.8f fac %s"%(str(zs[i]), pfloat(zs[i]), factor_sunit(zs[i])))

# ---- solve targets in the 8-dim quotient ----
xi = pickle.load(open('xi36_symbol.pkl','rb'))
t0v = Q.coord(xi['S0']); t1v = Q.coord(xi['S1'])
p0 = pi(t0v); p1 = pi(t1v)
sol0 = BM.gauss_jordan_solve(p0)[0]
sol1 = BM.gauss_jordan_solve(p1)[0]
assert BM*sol0 == p0 and BM*sol1 == p1
print("\nexact solution in quotient (coefficients c0=1-part, c1=sqrt3-part):")
for j,i in enumerate(basis_idx):
    print("   z=%-14s  c0=%-10s c1=%s"%(str(zs[i]), sol0[j], sol1[j]))
# lift back: verify full symbols match mod Decomp
comb0 = sp.zeros(27,1); comb1 = sp.zeros(27,1)
for j,i in enumerate(basis_idx):
    comb0 += sol0[j]*gens[i].T
    comb1 += sol1[j]*gens[i].T
diff0 = t0v.T - comb0; diff1 = t1v.T - comb1
# membership in Decomp:
rD = DM.rank()
print("diff0 in Decomp:", sp.Matrix.hstack(DM,diff0).rank()==rD, " diff1 in Decomp:", sp.Matrix.hstack(DM,diff1).rank()==rD)
pickle.dump({'basis_idx':basis_idx,'zs':[zs[i] for i in basis_idx],
             'sol0':[sp.Rational(x) for x in sol0], 'sol1':[sp.Rational(x) for x in sol1],
             'keep':keep}, open('xi36_li3_basis_solution.pkl','wb'))
print("saved xi36_li3_basis_solution.pkl")
