import sys, os, pickle, itertools
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

# ---- product symbols over V basis ('2','3','e') ----
def shuffle3(u, v, w):
    """symbol of product of three weight-1 logs with vecs u,v,w (multiplicity)."""
    S = {}
    for perm in set(itertools.permutations([(u,'u'),(v,'v'),(w,'w')])):
        pass
    # simpler: multiset-aware
    vecs = [u,v,w]
    for perm in itertools.permutations(range(3)):
        U,V,W = vecs[perm[0]], vecs[perm[1]], vecs[perm[2]]
        for a,ca in U.items():
            for b,cb in V.items():
                for c,cc in W.items():
                    S[(a,b,c)] = S.get((a,b,c),Fr(0)) + ca*cb*cc
    return {k:x for k,x in S.items() if x!=0}

def S_logprod(letters):
    """symbol of product log(L1)^e1 ... given list of basis-letter vecs with repetition"""
    vecs = [{l: Fr(1)} for l in letters]
    if len(vecs)==1:
        u = vecs[0]; return {(a,a,a): c**3 for a,c in u.items()} if False else tensor3(u,u,u)
    if len(vecs)==2:
        # log A * log B: A B + B A
        u,v = vecs
        S = {}
        for a,ca in u.items():
            for b,cb in v.items():
                S[(a,b,b)] = S.get((a,b,b),Fr(0))+ca*cb
                S[(b,a,b)] = S.get((b,a,b),Fr(0))+ca*cb
                S[(b,b,a)] = S.get((b,b,a),Fr(0))+ca*cb
        return {k:x for k,x in S.items() if x!=0}
    u,v,w = vecs
    return shuffle3(u,v,w)

def S_Li2log(w, p):
    """symbol of Li2(w)*log(p): w S-unit pair, p basis letter"""
    wv = vecfun(w); owv = vecfun(psub(ONE,w))
    S2 = {}
    for a,ca in owv.items():
        for b,cb in wv.items():
            S2[(a,b)] = S2.get((a,b),Fr(0)) - ca*cb
    pv = {p: Fr(1)}
    S = {}
    for (l1,l2_),c in S2.items():
        for t in [(p,l1,l2_),(l1,p,l2_),(l1,l2_,p)]:
            S[t] = S.get(t,Fr(0)) + c
    return {k:x for k,x in S.items() if x!=0}

# 19 visible products: 10 log cubes + 9 Li2*log
prod_syms = []; prod_names = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b
        letters = ['2']*a + ['3']*b + ['e']*c
        prod_syms.append(S_logprod(letters))
        prod_names.append("l2^%d l3^%d lu^%d"%(a,b,c))
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
for w in b2zs:
    for p in ('2','3','e'):
        prod_syms.append(S_Li2log(w,p))
        prod_names.append("Li2(%s)*l%s"%(str(w),p))
PM = sp.Matrix.hstack(*[ (lambda S: (lambda v: v)(sp.zeros(1,27))) (s) for s in []]) if False else None
Q = Quotient()
Pcols = [Q.coord(s) for s in prod_syms]
PM = sp.Matrix.hstack(*[c.T for c in Pcols])
print("rank of 19 product symbols:", PM.rank(), "(expect 19)")

# ---- self-test: T = Li21(x5, x6/x5) ----
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
T_sym = S_Li21(Xs[5], pdiv(Xs[6],Xs[5]))
# solve T_sym - sum a_i Li3 = products: use the earlier 8-basis solution
# recompute quickly the same 8 basis and sol as step_40
dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
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
solT = BM.gauss_jordan_solve(pi(Q.coord(T_sym)))[0]
# diff = S(T) - sum a Li3, solve over products
diff = Q.coord(T_sym).T
for j,i in enumerate(bidx):
    diff -= sp.Rational(solT[j])*gens[i].T
solP = PM.gauss_jordan_solve(diff)[0]
print("products exact solve ok:", PM*solP == diff)
if PM*solP != diff:
    # least: check membership
    print("diff in product span:", sp.Matrix.hstack(PM,diff).rank()==19)
else:
    print("exact product coefficients:")
    for j,nm in enumerate(prod_names):
        if solP[j]!=0: print("   %-28s %s"%(nm, sp.Rational(solP[j])))
pickle.dump({'prod_syms':prod_syms,'prod_names':prod_names,'bidx':bidx,'keep':keep,'indep':indep}, open('prod_sym_data.pkl','wb'))
