import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
from math import comb
import fastred
from fastred import *
from buildintegrand import build_integrand_pf
def build_tower(k):
    ipf=build_integrand_pf(k); work={}
    for key,R in ipf.items():
        kk,RR=fastred.fold_key(key,R); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    Hrec={}
    for deg in range(k,-1,-1):
        for key in sorted([kk for kk in work if sum(kk)==deg],reverse=True):
            R=work.pop(key)
            if deg==0: continue
            simple,higher,const=pf_split(R)
            Hrec[key]={r:e for r,e in enumerate(higher[1],start=1) if not e.is_zero()}
            if not higher[1] and not higher[2]: continue
            G=fastred.pf_integrate_higher(higher); C=pf_eval_at(G,BQ); Gm=[G[0],G[1],G[2],-C]
            i,j=key
            if i>=1:
                kk,RR=fastred.fold_key((i-1,j),pf_scale(pf_mul_DL(Gm,DL1),Fr(-i))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
            if j>=1:
                kk,RR=fold_key((i,j-1),pf_scale(pf_mul_DL(Gm,DL2),Fr(-j))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    return Hrec
def u(Hrec,key,delta):
    d=sum(key); return Hrec.get(key,{}).get(2*d-2-delta,ZERO)
def tau1(n):
    s=n+1; return ((-1)**s)*((IQ-BQ)**(-(s+1))+(IQ+AQ)**(-(s+1))-IQ**(-(s+1)))
def tau2(n):
    s=n+1; return ((-1)**s)*((IQ-AQ)**(-(s+1))+(IQ+BQ)**(-(s+1))-IQ**(-(s+1)))
def Fold_p(p,delta,Hrec):
    m=4*p-4-delta; tot=ZERO
    for mu in range(m,4*p-3):
        g=ZERO
        for dp in range(0,4*p-4-mu+1):
            g=g-u(Hrec,(p,p),dp)*tau1(4*p-4-mu-dp).scale(Fr(1,4*p-3-dp))
        tot=tot+g.conj_i()*(IQ**(-mu))*Fr(comb(mu-2,mu-m))*(IQ**(m-2))
    return tot.scale(-p)
def T1(p,q,k,delta,Hrec):
    d=p+q; tot=ZERO
    for dp in range(0,delta+1):
        tot = tot + u(Hrec,(p+1,q),dp)*tau1(delta-dp).scale(Fr(p+1,2*d-dp-1))
        if p>q:
            tot = tot + u(Hrec,(p,q+1),dp)*tau2(delta-dp).scale(Fr(q+1,2*d-dp-1))
    if q==p-1:
        tot = tot + Fold_p(p,delta,Hrec)
    return tot
bad=0; tested=0
for k in range(3,13):
    Hrec=build_tower(k)
    for d in range(2,k):
        for p in range(d,(d-1)//2,-1):
            q=d-p
            for delta in range(0,2*d-3):
                tested+=1
                if u(Hrec,(p,q),delta)!=T1(p,q,k,delta,Hrec):
                    bad+=1
                    if bad<=5: print("T1 MISMATCH",k,(p,q),delta)
print(f"FULL T1 recursion (my impl) vs data, k=3..12: {tested} cells, {bad} bad")
'''
with open(work+"step_17.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_17.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=900)
print(r.stdout); print("STDERR:",r.stderr[:2000])
with open(work+"progress.log","a") as fh: fh.write("step_17: FULL T1 recursion (my own impl) vs data, k=3..12\n")
