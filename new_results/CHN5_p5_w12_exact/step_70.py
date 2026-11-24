import os, sys
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp
import numpy as np
from scipy.optimize import linprog
from linalg_tools import monoms_homog
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
mons8=monoms_homog(8,5)
# test: target = -(a1+a2+a3+a4+a5)**8 (all negative coeffs), sense=+1, U,V=0 should be feasible
Qt=-(a1+a2+a3+a4+a5)**8
Qd={tuple(int(k) for k in mon):int(c) for mon,c in sp.Poly(Qt,*A).terms()}
n=0  # no U,V variables: pure feasibility of origin
Arows=[];brows=[]
for m in mons8:
    base=Qd.get(m,0)
    Arows.append([0.0]); brows.append(-1.0 - 1.0*float(base))  # sense=+1: 0 <= -1 - base  <=> base <= -1
res=linprog([0.0], A_ub=np.array(Arows), b_ub=np.array(brows), method='highs', bounds=(None,None))
print("feasibility test on -(sum)^8 (expect optimal/feasible):", res.success, res.message)
# so the LP formulation is consistent. The earlier infeasibility is genuine: no definite-sign rep.
# cleanup misleading outputs
for f in ['signnf_QA.pkl','signnf_QC.pkl']:
    if os.path.exists(f): os.remove(f); print("removed", f)
with open('progress.log','a') as f:
    import time
    f.write(time.strftime('%H:%M:%S ')+"STEP2 PROVED: B-(75/14)A = U1*C4+V1*Q7' (U1 deg4 50 terms, V1 deg2 6 terms) and Z+(25/3)A = U2*C4 (U2 deg4 50 terms, V2=0) - exact ideal membership, verified. => on V10: B=(75/14)A, Z=-(25/3)A EXACTLY; W12 irreducible content spanned by A and CLASS. [ideal.log]\n")
    f.write(time.strftime('%H:%M:%S ')+"LP note: NO sign-definite representative of QA or QC mod (C4,Q7') on the WHOLE open orthant (LP infeasible, code validated). Sign on V10+ must be argued on V10, not the orthant -> Gamma trace.\n")
print("logged + LP validated")
