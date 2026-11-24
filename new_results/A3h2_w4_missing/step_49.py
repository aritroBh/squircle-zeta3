import sys, pickle
from fractions import Fraction
from math import isqrt, gcd
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
r1=pickle.load(open(BASE+'/research/08_theory/A3h_weight4_target/rref_T3_1000000009.pkl','rb'))
r2=pickle.load(open(BASE+'/research/08_theory/A3h_weight4_target/rref_T3_1000000033.pkl','rb'))
p1,p2=1000000009,1000000033
M=p1*p2
def crt(a1,a2):
    return (a1+(a2-a1)*pow(p1,-1,p2)*p1)%M
def rat_recon(a,m):
    a%=m
    if a==0: return Fraction(0)
    r0,r1=m,a
    s0,s1=0,1
    bound=isqrt(m//2)
    while r1>bound:
        q=r0//r1
        r0,r1=r1,r0-q*r1
        s0,s1=s1,s0-q*s1
    if s1<0: r1,s1=-r1,-s1
    if gcd(r1,s1)!=1 or s1==0: return None
    return Fraction(r1,s1)
# test entries across the matrix
import random
random.seed(0)
tested=0; ok=0; vals=set(); fails=0
for i in range(0,5074,250):
    row1=dict(r1['rowsF'][i]); row2=dict(r2['rowsF'][i])
    for j,v1 in list(row1.items()):
        v2=row2.get(j,0)
        x=crt(v1,v2)
        val=rat_recon(x,M)
        tested+=1
        if val is None: fails+=1
        else:
            ok+=1; vals.add(val)
print(f'tested {tested}, reconstructed {ok}, failed {fails}, distinct values {len(vals)}')
print('distinct values (first 30):', sorted(vals, key=lambda f: (f.denominator,f.numerator))[:30])
print('max denominator:', max((v.denominator for v in vals), default=None))
