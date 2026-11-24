import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
Q = Quotient()

# 19 independent Decomp columns
DM_full = sp.Matrix.hstack(*[d.T for d in Q.D])
indep = []
M = sp.zeros(27,0)
for j in range(DM_full.cols):
    M2 = sp.Matrix.hstack(M, DM_full[:,j]) if M.cols else DM_full[:,j]
    if M.cols==0 or M2.rank()>M.rank():
        indep.append(j); M = M2
    if M.rank()==19: break
D19 = M
print("indep decomp cols:", len(indep), "rank", D19.rank())
keep = []
for i in range(27):
    e = sp.zeros(27,1); e[i]=1
    M2 = sp.Matrix.hstack(M, e)
    if M2.rank() > M.rank():
        keep.append(i); M = M2
    if len(keep)==8: break
assert M.rank()==27
A = M   # 27x27: [D19 | E8]
Ainv = A.inv()
def pi(v):
    w = Ainv * v.T
    return sp.Matrix(w[19:,:])
Gp = [pi(g) for g in gens]
def niceness(i):
    z = zs[i]; _,v2,v3,k = factor_sunit(z)
    fz = pfloat(z)
    pref = 0 if 0 < fz < 1 else (1 if fz < 0 else 2)
    return (pref, abs(v2)+abs(v3)+abs(k))
order = sorted(range(len(zs)), key=niceness)
basis_idx = []
BM = sp.zeros(8,0)
for i in order:
    M2 = sp.Matrix.hstack(BM, Gp[i]) if BM.cols else Gp[i]
    if BM.cols==0 or M2.rank()>BM.rank():
        basis_idx.append(i); BM = M2
    if len(basis_idx)==8: break
print("basis gens:", basis_idx, "rank:", BM.rank())
for i in basis_idx:
    print("   z = %-14s ~ %+.8f fac %s"%(str(zs[i]), pfloat(zs[i]), factor_sunit(zs[i])))

xi = pickle.load(open('xi36_symbol.pkl','rb'))
t0v = Q.coord(xi['S0']); t1v = Q.coord(xi['S1'])
p0 = pi(t0v); p1 = pi(t1v)
sol0 = BM.gauss_jordan_solve(p0)[0]
sol1 = BM.gauss_jordan_solve(p1)[0]
print("\nexact solution (c0 = 1-part, c1 = sqrt3-part):")
tot0 = sp.zeros(27,1); tot1 = sp.zeros(27,1)
for j,i in enumerate(basis_idx):
    print("   z=%-14s  c0=%-10s c1=%s"%(str(zs[i]), sp.Rational(sol0[j]), sp.Rational(sol1[j])))
    tot0 += sp.Rational(sol0[j])*gens[i].T
    tot1 += sp.Rational(sol1[j])*gens[i].T
diff0 = t0v.T - tot0; diff1 = t1v.T - tot1
print("diff0 in Decomp:", sp.Matrix.hstack(D19,diff0).rank()==19,
      " diff1 in Decomp:", sp.Matrix.hstack(D19,diff1).rank()==19)
pickle.dump({'basis_idx':basis_idx,'zs':[zs[i] for i in basis_idx],
             'sol0':[sp.Rational(x) for x in sol0], 'sol1':[sp.Rational(x) for x in sol1],
             'keep':keep,'indepD':indep}, open('xi36_li3_basis_solution.pkl','wb'))
print("saved xi36_li3_basis_solution.pkl")
