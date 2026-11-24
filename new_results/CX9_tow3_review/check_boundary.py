"""Exact first-level, E3 and lowest admissible degree checks against TOWN1."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
from functools import lru_cache
import pickle,sys
import sympy as sp
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'TOW3_tower_depth6'))
sys.path.insert(0,str(ROOT/'WTS_weight_k_structure'/'theory'))
import tow3lib as T
from fastred import QSI,ZERO,IQ,AQ,BQ
H=pickle.load(open(ROOT/'TOWN1_tower_large_n'/'data'/'towers_D8_K52.pkl','rb'))
def qsi(expr):
    z=sp.expand(expr)
    return QSI(F(sp.re(z)),0,F(sp.im(z)),0)
@lru_cache(None)
def tau(n,which):
    s=n+1
    return (-1)**s*((IQ-(BQ if which==1 else AQ))**(-s-1)+(IQ+(AQ if which==1 else BQ))**(-s-1)-IQ**(-s-1))
for n in range(8):
    assert tau(n,1)==qsi(T.tau1(n)) and tau(n,2)==qsi(T.tau2(n))
print('PASS tau n=0..7 from DL definitions equals submitted table',flush=True)
def actual(key,k,dp): return QSI(*H[key][dp][k])
@lru_cache(None)
def off(key,dp,s):
    return actual(key,sum(key),0)*qsi(T.e(dp,s,sum(key),key[0]))
@lru_cache(None)
def hat(p,dp,s):
    return sum((off((p+1,p),j,s)*tau(dp-j,1)*F(p+1,4*p-j-1) for j in range(dp+1)),ZERO)
def stable(p,k,dp):
    # Generalized binomial at -1 is the actual Newton extrapolation to k=2p.
    n=k-2*p-1
    return sum((hat(p,dp,s)*((-1)**s if n==-1 else comb(n,s)) for s in range(dp//2+1)),ZERO)
def candidate(key,k,dp):
    p,q=key; d=p+q
    if p==q: return off(key,dp,0) if k==d else stable(p,k,dp)
    return sum((off(key,dp,s)*comb(k-d,s) for s in range(dp//2+1)),ZERO)
def diagfeed(p,delta,U):
    direct=sum((U(j)*tau(delta-j,2)*F(p,4*p-3-j) for j in range(delta+1)),ZERO)
    ft=ZERO; r=4*p-4-delta
    for mu in range(r,4*p-3):
        g=-sum((U(j)*tau(4*p-4-mu-j,1)*F(1,4*p-3-j) for j in range(4*p-4-mu+1)),ZERO)
        ft+=g.conj_i()*IQ**(-mu)*comb(mu-2,mu-r)*IQ**(r-2)*(-p)
    return direct+ft
count=0; discrepancy=0
for p in (3,4,5,6):
    for key in ((p,p-1),(p,p)):
        d=sum(key)
        for k in range(d,d+6):
            for dp in range(7):
                assert candidate(key,k,dp)==actual(key,k,dp),(key,k,dp)
                count+=1
    for delta in (5,6):
        diff=lambda dp:actual((p,p),2*p,dp)-stable(p,2*p,dp)
        assert diagfeed(p,delta,diff).is_zero(),(p,delta)
        discrepancy+=not diff(delta).is_zero()
print(f'PASS {count} exact boundary tower cells, p=3..6, diagonal/adjacent, first 6 levels, depths 0..6',flush=True)
print(f'PASS 8 E3 direct+fold checks at p=3..6, delta=5,6; {discrepancy}/8 diagonal initial discrepancies nonzero',flush=True)
# d=5 is indispensable for weight five, including delta=6.
count=0
for key in ((5,0),(4,1),(3,2)):
    for k in (5,6,7,8,9,12,21,40,52):
        for dp in (5,6):
            assert candidate(key,k,dp)==actual(key,k,dp),(key,k,dp)
            count+=1
print(f'PASS {count} d=5 depth-5/6 cells, all three keys, levels 5,6,7,8,9,12,21,40,52',flush=True)
p=sp.Symbol('p',integer=True,positive=True)
for delta in (5,6):
    s=delta//2
    top=sp.factor((p+1)*sum(T.tau1(delta-j)*T.e(j,s,2*p+1,p+1)*2*sp.I*(4*p-1)/(p+1)/(4*p-j-1) for j in range(delta+1)))
    print(f'normalized diagonal top ehat({delta},{s})/rho_pp(0) =',top,flush=True)
