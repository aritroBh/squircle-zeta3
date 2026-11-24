import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g=S.init()
d=pickle.load(open('delta14.pkl','rb'))
sol_c=d['sol_c']; CRs=d['CRs']
LI3=[S.li3_symbol(z) for z in CRs]
letters=[0,1,2]; basis3=list(itertools.product(letters,repeat=3)); bidx={t:i for i,t in enumerate(basis3)}
Dg=[]
for a in letters:
    for b in letters:
        for c in letters:
            v={}
            for t in [(a,b,c),(b,a,c),(b,c,a)]: v[t]=v.get(t,0)+1
            Dg.append(v)
def t2v(T):
    v=[Fr(0)]*27
    for t,c in T.items():
        if t in bidx: v[bidx[t]]+=Fr(c)
    return v
Dgv=[t2v(T) for T in Dg]; LI3v=[t2v(T) for T in LI3]
Mcols=LI3v+Dgv
M=sp.Matrix(27,len(Mcols),lambda i,j: sp.Rational(Mcols[j][i]))
sigma=pickle.load(open('sigma14.pkl','rb'))
# for each of 3 components: solve fully (keep all coeffs), check residual exact 0
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,dd=cFT
    if key in bidx:
        s_comp[0][bidx[key]]+=Fr(A)/dd; s_comp[1][bidx[key]]+=Fr(B)/dd; s_comp[2][bidx[key]]+=Fr(C)/dd
for k in range(3):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    res=M.gauss_jordan_solve(tgt); sol=res[0]
    free=set()
    for z in range(sol.rows): free|=sol[z].free_symbols
    sol=sol.subs({fs:0 for fs in free})
    resid = M*sol - tgt
    print(f"component {k}: exact residual = {sp.simplify(resid).is_zero_matrix if hasattr(resid,'is_zero_matrix') else all(x==0 for x in resid)}")
# stuffle check for S_Li21 (Lam convention): S(Li21(x,y))+S(Li12(y,x))+S(Li3(xy)) = S(Li2(x))*S(Li1(y)) (symbol)
# use x=x1, y=x2/x1. Li12(y,x) = swap. Compute via sym_Lam formulas.
def sym_Li21(x,y):
    L1=S.vec(E.fsub(E.fmul(x,y),E.ONE)); L2=S.vec(E.fmul(x,y)); L3=S.vec(E.fsub(y,E.ONE))
    L4=S.vec(E.fdiv(E.fmul(x,E.fsub(y,E.ONE)),E.fsub(x,E.ONE))); L5=S.vec(E.fsub(x,E.ONE)); L6=S.vec(x)
    Sdod=S.sadd(S.tensor(L1,L2,L3),S.tensor(L1,L4,L6),S.tensor(L3,L5,L6))
    out={}
    for (a,b,c),cf in Sdod.items():
        out[(c,b,a)]=out.get((c,b,a),Fr(0))-cf
        if out[(c,b,a)]==0: del out[(c,b,a)]
    return out
def sym_Li12(x,y):
    # Li_{1,2}(x,y): stuffle partner. Its symbol = swap of Li21? Use S(Li12(x,y)) from Li21(y,x)? 
    # Li_{1,2}(x,y) = prefix f(y,x). Compute via S_Li21 swapped formula: Li_{1,2}(x,y) letters...
    # standard: S(Li_{1,2}(x,y)) = tau of S(Li_{2,1}(y,x))? Just compute via DOD formula for Li_{1,2}.
    pass
# Simpler convention check: verify sym_Li21 matches the S3 recursion on the G-word I(0;1/x_j,0,1/x_i;1)
x1=E.X[1]; x2=E.X[2]
lam12_sym21 = sym_Li21(x1, E.fdiv(x2,x1))
lam12_s3 = S.S3(E.finv(x2), E.ZERO, E.finv(x1))   # I(0; 1/x2, 0, 1/x1; 1)
# compare (up to overall sign convention?) 
def ssub(A,B):
    R=dict(A)
    for k,c in B.items():
        R[k]=R.get(k,Fr(0))-c
        if R[k]==0: del R[k]
    return R
d1=ssub(lam12_sym21, lam12_s3)
print("Lam(1,2): S_Li21 vs S3(1/x2,0,1/x1) difference terms:", len(d1))
if len(d1)<15: print("  diff:", d1)
