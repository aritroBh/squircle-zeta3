
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
from math import comb
import fastred
from fastred import *
from buildintegrand import build_integrand_pf

# build full work-items W (simple+higher) AND Gm for each key at level k, and tower u
def reduce_full(k):
    ipf=build_integrand_pf(k); work={}
    for key,R in ipf.items():
        kk,RR=fastred.fold_key(key,R); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    Gm={}; Hrec={}
    for deg in range(k,-1,-1):
        for key in sorted([kk for kk in work if sum(kk)==deg],reverse=True):
            R=work.pop(key)
            if deg==0: continue
            simple,higher,const=pf_split(R)
            Hrec[key]={r:e for r,e in enumerate(higher[1],start=1) if not e.is_zero()}
            if not higher[1] and not higher[2]: continue
            G=fastred.pf_integrate_higher(higher); C=pf_eval_at(G,BQ); Gm[key]=[G[0],G[1],G[2],-C]
            i,j=key
            if i>=1:
                kk,RR=fastred.fold_key((i-1,j),pf_scale(pf_mul_DL(Gm[key],DL1),Fr(-i))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
            if j>=1:
                kk,RR=fold_key((i,j-1),pf_scale(pf_mul_DL(Gm[key],DL2),Fr(-j))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    return Gm,Hrec
def u_of(Hrec,key,delta):
    d=sum(key); return Hrec.get(key,{}).get(2*d-2-delta,ZERO)

# T1 fold formula (my own): Fold_p(delta; u_pp), u_pp from Hrec
def Fold_p(p,delta,Hrec,k):
    m=4*p-4-delta; tot=ZERO
    for mu in range(m,4*p-3):
        g=ZERO
        for dp in range(0,4*p-4-mu+1):
            g = g - u_of(Hrec,(p,p),dp)*tau1(4*p-4-mu-dp).scale(Fr(1,4*p-3-dp))
        tot = tot + g.conj_i()*(IQ**(-mu))*Fr(comb(mu-2,mu-m))*(IQ**(m-2))
    return tot.scale(-p)
def tau1(n):
    s=n+1; return ((-1)**s)*((IQ-BQ)**(-(s+1))+(IQ+AQ)**(-(s+1))-IQ**(-(s+1)))

# actual fold feed: +i order-(2d-2-delta) coeff of -p*F(Gm_pp DL1), d=2p-1
def actual_fold(p,delta,Gm,k):
    d=2*p-1
    FF = pf_fold(pf_scale(pf_mul_DL(Gm[(p,p)],DL1),Fr(-p)))   # -p F(Gm DL1)
    # +i principal part coefficient at order (2d-2-delta)
    ipos = FF[1]; r = 2*d-2-delta
    return ipos[r-1] if r-1 < len(ipos) else ZERO

bad=0; tested=0
for k in range(4,13):
    Gm,Hrec = reduce_full(k)
    for p in range(2, (k//2)+1):
        if 2*p-1 > k: continue
        if (p,p) not in Gm: continue
        for delta in range(0, min(4, 2*(2*p-1)-4)+1):
            tested+=1
            mine = Fold_p(p,delta,Hrec,k)
            act = actual_fold(p,delta,Gm,k)
            if mine!=act:
                bad+=1
                if bad<=5: print("FOLD MISMATCH",k,p,delta, mine.sympy(), act.sympy())
print(f"fold formula Fold_p vs actual reducer fold, k=4..12: {tested} cells, {bad} bad")
