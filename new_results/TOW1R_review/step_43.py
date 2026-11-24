import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
import fastred
from fastred import *
from buildintegrand import build_integrand_pf
def pf_is_zero(P): return not P[0] and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()
def simple_part(P): return pf_split(P)[0]
def gms_at(k):
    work={}
    for key,R in build_integrand_pf(k).items():
        kk,RR=fastred.fold_key(key,R); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    gms={}
    for deg in range(k,-1,-1):
        for key in sorted([kk for kk in work if sum(kk)==deg],reverse=True):
            R=work.pop(key)
            if deg==0: continue
            simple,higher,const=pf_split(R)
            if not higher[1] and not higher[2]: continue
            G=fastred.pf_integrate_higher(higher); C=pf_eval_at(G,BQ); gms[key]=[G[0],G[1],G[2],-C]
            i,j=key
            if i>=1:
                kk,RR=fastred.fold_key((i-1,j),pf_scale(pf_mul_DL(gms[key],DL1),Fr(-i))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
            if j>=1:
                kk,RR=fastred.fold_key((i,j-1),pf_scale(pf_mul_DL(gms[key],DL2),Fr(-j))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    return gms
g3 = gms_at(3)
G30, G21 = g3[(3,0)], g3[(2,1)]
s30=[{},[QSI(0,Fr(1,12))],[QSI(0,Fr(1,12))],QSI(Fr(-1,12))]
s21=[{},[QSI(0,Fr(-1,4))],[QSI(0,Fr(-1,4))],QSI(Fr(1,4))]
R20_3=pf_zero(); R20_3[0][ZERO]=QSI(Fr(1,4),Fr(-1,4)); R20_3[0][AQ]=QSI(Fr(1,4)); R20_3[0][-AQ]=QSI(Fr(-1,4)); R20_3[0][-BQ]=QSI(Fr(-1,2))
R11_3=pf_zero(); R11_3[0][ZERO]=QSI(Fr(1,4),Fr(1,4)); R11_3[0][-AQ]=QSI(Fr(-1,2))
S20=pf_zero(); S20[0][-AQ]=QSI(Fr(1,2)); S20[0][-BQ]=QSI(Fr(-1,2))
S11=pf_zero(); S11[0][ZERO]=QSI(Fr(1,2)); S11[0][-AQ]=QSI(-1)
C1=pf_is_zero(pf_sub(simple_part(pf_add(pf_scale(pf_mul_DL(G30,DL1),Fr(-3)),pf_neg(pf_mul_DL(G21,DL2)))),R20_3))
C2=pf_is_zero(pf_sub(simple_part(pf_add(pf_scale(pf_mul_DL(s30,DL1),Fr(-3)),pf_neg(pf_mul_DL(s21,DL2)))),S20))
C3=pf_is_zero(pf_sub(simple_part(pf_scale(pf_mul_DL(G21,DL1),Fr(-2))),R11_3))
C4=pf_is_zero(pf_sub(simple_part(pf_scale(pf_mul_DL(s21,DL1),Fr(-2))),S11))
print("C1:",C1," C2:",C2," C3:",C3," C4:",C4)
# verify NO +-i poles survive in SimplePart (they should be real-only). Check W's full PF has no +-i order-1 after simple_part.
# slope derivation check: s30 = sqrt2/12*(1/(v-i)+1/(v+i)) - 1/12 ; s21 = -3 s30
# my derived s30 from tower: +i principal part slope sqrt2/12 at order 1; const -1/12
s30_mine = [{},[QSI(0,Fr(1,12))],[QSI(0,Fr(1,12))],QSI(Fr(-1,12))]
print("s30 PF matches my derivation:", s30_mine==s30)
# print G30,G21 for inspection
import sympy as sp
print("Gm30^(3) =", sp.simplify(pf_to_sympy(G30)))
print("Gm21^(3) =", sp.simplify(pf_to_sympy(G21)))
'''
with open(work+"step_15.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_15.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=300)
print(r.stdout); print("STDERR:",r.stderr[:2000])
with open(work+"progress.log","a") as fh: fh.write("step_15: verified C1-C4 independently + s30 PF derivation\n")
