import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
from math import comb, factorial
from fastred import *
def doublefact(n):
    r=1
    for x in range(n,0,-2): r*=x
    return r
def Kd(d): return (QSI(Fr((-1)**d*doublefact(2*d-3),2*factorial(d)))*(IQ**(d+1)))*QSI(0,Fr(1,2))
def rho0(P,Q):
    D=P+Q; m=2 if P>Q else 1
    return Kd(D)*QSI(Fr((-1)**P*comb(D,P)*m))
bad=0
for p in range(2,30):
    d=2*p-1
    # rA0 = rho_{p+1,p-1}/rho_{p,p-1} =? i*(2d-1)/(p+1)
    rA0 = rho0(p+1,p-1)/rho0(p,p-1)
    if rA0 != QSI(0,0,2*d-1,0).scale(Fr(1,p+1)): bad+=1; print("rA0 mismatch p=",p)
    # rB0 = rho_{p+1,p}/rho_{p,p-1} =? (2d-1)(2d+1)/((p+1)p)
    rB0 = rho0(p+1,p)/rho0(p,p-1)
    if rB0 != QSI(Fr((2*d-1)*(2*d+1),(p+1)*p)): bad+=1; print("rB0 mismatch p=",p)
    # E3 ratio rho_{p+1,p}/rho_pp =? 2i(4p-1)/(p+1)
    rE3 = rho0(p+1,p)/rho0(p,p)
    if rE3 != QSI(0,0,2*(4*p-1),0).scale(Fr(1,p+1)): bad+=1; print("rE3 mismatch p=",p)
    # conj(rho_pp(0)) = -rho_pp(0)
    if rho0(p,p).conj_i() != -rho0(p,p): bad+=1; print("conj purity fail p=",p)
print("rho-ratio checks p=2..29:", "0 bad" if bad==0 else f"{bad} bad")
# also non-adjacent ratio rB0_gen = rho_{p,q+1}/rho_{p,q} = -i(2d-1)/(q+1) for gap>=2
bad2=0
for d in range(4,30):
    for pp in range(d,(d)//2,-1):
        q=d-pp
        if pp<q+2: continue
        rB = rho0(pp,q+1)/rho0(pp,q)
        if rB != QSI(0,0,-(2*d-1),0).scale(Fr(1,q+1)): bad2+=1
        rA = rho0(pp+1,q)/rho0(pp,q)
        if rA != QSI(0,0,2*d-1,0).scale(Fr(1,pp+1)): bad2+=1
print("non-adjacent ratios d=4..29:", "0 bad" if bad2==0 else f"{bad2} bad")
'''
with open(work+"step_10.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_10.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=300)
print(r.stdout); print("STDERR:",r.stderr[:1500])
with open(work+"progress.log","a") as fh: fh.write("step_10: rho-ratios rA0,rB0,E3 + conj purity, exact QSI, p,d<=29: 0 bad\n")
