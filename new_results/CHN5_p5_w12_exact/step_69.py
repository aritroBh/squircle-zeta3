import os, sys, time, pickle
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
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
mons8=monoms_homog(8,5)
mons4=monoms_homog(4,5); mons2=monoms_homog(2,5)
C4d=sym_to_dict(C4); Q7d=sym_to_dict(Q7)
nU=len(mons4); nV=len(mons2); n=nU+nV
def feas(Qtarget, sense):
    """find U,V: sense*(Q - U.C4 - V.Q7)_m <= -1 for all m (float), then exact-verify."""
    Qd=sym_to_dict(Qtarget)
    Arows=[]; brows=[]
    for m in mons8:
        row=[0.0]*n
        base=Qd.get(m,F(0))
        for i,um in enumerate(mons4):
            rem=tuple(m[k]-um[k] for k in range(5))
            if all(r>=0 for r in rem) and rem in C4d:
                row[i] = -sense*float(C4d[rem])
        for j,vm in enumerate(mons2):
            rem=tuple(m[k]-vm[k] for k in range(5))
            if all(r>=0 for r in rem) and rem in Q7d:
                row[nU+j] = -sense*float(Q7d[rem])
        Arows.append(row); brows.append(-1.0 - sense*float(base))
    res=linprog([0.0]*n, A_ub=np.array(Arows), b_ub=np.array(brows), method='highs', bounds=(None,None))
    return res, Qd
def exact_round_verify(Qtarget, xsol, sense, scale=10**7):
    Ur=[int(round(xsol[i]*scale)) for i in range(nU)]
    Vr=[int(round(xsol[nU+j]*scale)) for j in range(nV)]
    U=sum(sp.Integer(Ur[i])*sp.prod(A[k]**mons4[i][k] for k in range(5)) for i in range(nU) if Ur[i])
    V=sum(sp.Integer(Vr[j])*sp.prod(A[k]**mons2[j][k] for k in range(5)) for j in range(nV) if Vr[j])
    R=sp.expand(Qtarget*sp.Integer(scale)**1 - 0) # placeholder
    # R_exact = scale*Qtarget - U*C4 - V*Q7  (so R_exact/scale ~ R_float)
    R=sp.expand(sp.Integer(scale)*Qtarget - sp.expand(U*C4) - sp.expand(V*Q7))
    Rd=sym_to_dict(R)
    coeffs=[Rd.get(m,F(0)) for m in mons8]
    badpos=sum(1 for c in coeffs if sense*c>=0)
    nzero=sum(1 for c in coeffs if c==0)
    minc=min((sense*c for c in coeffs), default=0)
    return R, Rd, badpos, nzero, minc
with open('sign_nf2.log','w') as lg:
    for name,Qt,sense in [('QA<0 (=> A<0 on V10+)',QA,+1), ('QC>0 (=> CLASS>0 on V10+)',QC,-1), ('QC<0',QC,+1)]:
        res,Qd=feas(Qt,sense)
        if not res.success:
            lg.write(f"{name}: LP infeasible/failed: {res.message}\n"); lg.flush(); continue
        R,Rd,badpos,nzero,minc=exact_round_verify(Qt,res.x,sense)
        lg.write(f"{name}: float feas OK. Exact rounded: coeffs violating strict sign={badpos}, zeros={nzero}, min sense*coeff={minc}\n")
        if badpos==0 and nzero<495:
            pickle.dump({'R':str(R),'sense':sense}, open(f'signnf_{"QA" if name.startswith("QA") else "QC"}.pkl','wb'))
            lg.write(f"{name}: *** CERTIFIED definite-sign representative (deg8, all {495-nzero} nonzero coeffs strict sign) ***\n")
        lg.flush()
    lg.write("DONE\n")
print(open(W+'/sign_nf2.log').read())
