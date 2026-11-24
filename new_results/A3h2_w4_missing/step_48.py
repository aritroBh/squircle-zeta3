import sys, pickle, json
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
from fractions import Fraction as F
# check fw ordering
fwj=json.load(open(BASE+'/research/08_theory/A3h_weight4_target/T3_quotient_nf.json'))['freewords']
d=pickle.load(open(WORK+'/cands29_reduced.pickle','rb'))
print('fw ordering matches json:', [list(w) for w in d['freewords']]==fwj)
# words(4) to map free cols
sys.path.insert(0,BASE+'/research/08_theory/A3h_weight4_target/Q')
from relations import words
cols=words(4)
r1=pickle.load(open(BASE+'/research/08_theory/A3h_weight4_target/rref_T3_1000000009.pkl','rb'))
free_words=[cols[j] for j in r1['free']]
print('free words match:', [list(w) for w in free_words]==fwj)
# CRT + rational reconstruction test on rowsF entries
r2=pickle.load(open(BASE+'/research/08_theory/A3h_weight4_target/rref_T3_1000000033.pkl','rb'))
p1,p2=1000000009,1000000033
assert r1['free']==r2['free'] and len(r1['rowsF'])==len(r2['rowsF'])
import sympy as sp
def crt(a1,p1,a2,p2):
    m=p1*p2
    return (a1+(a2-a1)*pow(p1,-1,p2)*p1)%m
M=p1*p2
def ratrec(x,M):
    try:
        return sp.nsimplify(x,(M,),rational=True) if False else sp.Rational(sp.symmetric_residue(x,M)[0] if False else x)
    except Exception: return None
from sympy.polys.domains import ZZ
def rational_reconstruction(a,m):
    # find p/q with p/q = a mod m, |p|,q small
    from sympy import Rational
    r0,r1=m,a%m
    s0,s1=0,1
    bound=sp.isqrt(m)//2
    while r1:
        q=(r0//r1)
        r0,r1=r1,r0-q*r1
        s0,s1=s1,s0-q*s1
        if r1 and r1*r1<=m and abs(s1)*s1<=m//2:
            pass
    # use sympy
    return None
# use sympy's nsimplify-free method: fractions via sympy.polys
from sympy import floor, Rational
def rr(a,m):
    """rational reconstruction a mod m -> p/q with |p|,q <= sqrt(m/2)"""
    a%=m
    if a==0: return Rational(0)
    u0,u1,u2=1,0,m
    v0,v1,v2=0,1,a
    while v2:
        q=u2//v2
        u0,u1,u2,v0,v1,v2=v0,v1,v2,u0-q*v0,u1-q*v1,u2-q*v2
        if u2 and abs(v2) and v2*v2<=m and v1*v1*2<=m:
            pass
    # fallback simple loop
    return None
# simplest: sympy.simplify via continued fractions using sympy.ntheory
from sympy.ntheory.modular import solve_congruence
def rat_recon(a,m):
    a%=m
    if a==0: return Rational(0,1)
    x,y=m,a
    A,B=0,1
    while y:
        q=x//y
        x,y=y,x-q*y
        A,B=B,A-q*B
        if x and B and x*x<=m and 2*B*B<=m:
            return Rational((-1)**0 * x * (1 if (a*B)%m==x%m else -1), abs(B))
    return None
# test on first 20 entries of rowsF
import time
t=time.time()
count=0; big=0
samples=[]
for i in list(range(0,5074,500)):
    for j,v1 in r1['rowsF'][i][:5]:
        # find matching entry in r2
        d2=dict(r2['rowsF'][i])
        v2=d2.get(j,0)
        x=crt(v1,p1,v2,p2)
        val=rat_recon(x,M)
        if val is None:
            big+=1
        else:
            samples.append(val)
        count+=1
print('tested',count,'entries; failed recon:',big,'time',time.time()-t)
print('sample reconstructed values:',samples[:15])
