import os
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
signnf = r'''
import sys, time, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
import numpy as np
from scipy.optimize import linprog
from linalg_tools import monoms_homog
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
def dict_to_sym(p):
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator*sp.prod(A[i]**e[i] for i in range(5)) for e,c in p.items()))
def sym_to_dict(poly):
    P=sp.Poly(sp.expand(poly),*A)
    return {tuple(int(k) for k in mon): F(int(c.p),int(c.q)) for mon,c in P.terms()}
Q=pickle.load(open('Q8_quotients.pkl','rb'))
QA=dict_to_sym(Q['A']); QC=dict_to_sym(Q['CLASS'])
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
mons8=monoms_homog(8,5); midx={m:i for i,m in enumerate(mons8)}
mons4=monoms_homog(4,5); mons2=monoms_homog(2,5)
C4d=sym_to_dict(C4); Q7d=sym_to_dict(Q7)
nU=len(mons4); nV=len(mons2); n=nU+nV
# columns: variable = coeff of U at mons4[i], V at mons2[j]. R_m = Q_m - sum U_i C4_{m-i} - sum V_j Q7_{m-j}
def build(Qtarget, sense):
    """minimize t s.t. for all m: sense*(Q_m - U.C4 - V.Q7)_m <= t. sense=+1: want R_m<0; -1: want R_m>0."""
    Qd=sym_to_dict(Qtarget)
    # variables x = [U_0..U_{nU-1}, V_0..V_{nV-1}, t]
    Arows=[]; brows=[]
    for m in mons8:
        row=[0.0]*(n+1)
        base=Qd.get(m,F(0))
        for i,um in enumerate(mons4):
            rem=tuple(m[k]-um[k] for k in range(5))
            if all(r>=0 for r in rem) and rem in C4d:
                row[i] = -sense*float(C4d[rem])
        for j,vm in enumerate(mons2):
            rem=tuple(m[k]-vm[k] for k in range(5))
            if all(r>=0 for r in rem) and rem in Q7d:
                row[nU+j] = -sense*float(Q7d[rem])
        row[n]=1.0
        Arows.append(row); brows.append(sense*float(base))
    c=[0.0]*n+[1.0]
    res=linprog(c, A_ub=np.array(Arows), b_ub=np.array(brows), method='highs', bounds=[(None,None)]*n+[(None,None)])
    return res, Qd
def verify(Qtarget, sol, sense):
    """exact: R = Qtarget - U*C4 - V*Q7'; check sense*R has all coeff < 0 (and none zero)."""
    U=sum(sp.Rational(sol[i]).limit_denominator(10**9)*sp.prod(A[k]**mons4[i][k] for k in range(5)) for i in range(nU) if abs(sol[i])>1e-12)
    V=sum(sp.Rational(sol[nU+j]).limit_denominator(10**9)*sp.prod(A[k]**mons2[j][k] for k in range(5)) for j in range(nV) if abs(sol[nU+j])>1e-12)
    R=sp.expand(Qtarget - sp.expand(U*C4) - sp.expand(V*Q7))
    Rd=sym_to_dict(R)
    coeffs=[Rd.get(m,F(0)) for m in mons8]
    nbad=sum(1 for c in coeffs if sense*c>=0)
    nzero=sum(1 for c in coeffs if c==0)
    return R, Rd, nbad, nzero, len([c for c in coeffs if c!=0])
with open('sign_nf.log','w') as lg:
    for name,Qt,sense in [('QA (want R<0 => A<0)',QA,+1), ('QC (want R>0 => CLASS>0)',QC,-1), ('QC (want R<0 => CLASS<0)',QC,+1)]:
        res,Qd=build(Qt,sense)
        if not res.success:
            lg.write(f"{name}: LP failed: {res.message}\n"); lg.flush(); continue
        t=res.fun
        lg.write(f"{name}: LP min max-coeff t = {t} (float)\n"); lg.flush()
        if t<0:
            R,Rd,nbad,nzero,nnz=verify(Qt,res.x,sense)
            lg.write(f"{name}: exact verify: coeffs violating strict sign = {nbad}, zeros = {nzero}, nonzero coeffs = {nnz}/495\n")
            if nbad==0:
                pickle.dump({'R':str(R)}, open(f'signnf_{"QA" if "QA" in name else "QC"}.pkl','wb'))
                lg.write(f"{name}: *** CERTIFIED definite-sign normal form ***\n")
        lg.flush()
    lg.write("DONE\n")
print("sign nf done")
'''
open(W+'/run_sign_nf.py','w').write(signnf)
print("run_sign_nf.py written")
