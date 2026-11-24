import mpmath as mp, json, os
from mpmath import mpf
mp.mp.dps=65
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
print("bgA.log tail:")
print(open(WORK+'/bgA.log').read()[-400:])
z3=mp.zeta(3); pi=mp.pi; l2=mp.log(2); l3=mp.log(3); sq5=mp.sqrt(5)
def sermul(a,b,K):
    o=[mpf(0)]*(K+1)
    for i in range(K+1):
        for j in range(i+1): o[i]+=a[j]*b[i-j]
    return o
def serinv(a,K):
    o=[mpf(0)]*(K+1); o[0]=1/a[0]
    for n in range(1,K+1): o[n]=-sum(a[j]*o[n-j] for j in range(1,n+1))/a[0]
    return o
L=(mpf(-6)/5)*l2+(mpf(13)/5)*l3
u=[mpf(0),1,-1,1,-1];Lu=[L*x for x in u];e=[mpf(0)]*5;e[0]=1;pw=[mpf(0)]*5;pw[0]=1
for k in range(1,5):
    pw=sermul(pw,Lu,4)
    for i in range(4,-1,-1): e[i]+=pw[i]/mp.factorial(k)
finv=serinv(e,4)
SA=[mpf(s) for s in json.load(open(WORK+'/modelA_h128.json'))['Sdk']]
T=sermul(list(SA),finv,4)
t2n=T[2]/sq5
def try_pslq(tag,target, basis, tol=mpf('1e-45')):
    names=list(basis); x=[target]+[basis[n] for n in names]
    rel,resid=pslq(x,tol=tol,maxsteps=700)
    if rel is None: print("   %s: none"%tag); return
    h=max(abs(r) for r in rel)
    print("   %s: h=%.2e resid=%.1e rel=%s"%(tag,float(h),float(resid),rel))
# Hypothesis b: weight-2 = {1, pi2, ln^2, ln ln} only
try_pslq("{1,pi2,ln2^2,ln3^2,ln2ln3}", t2n, {'1':mpf(1),'pi2':pi**2,'l2^2':l2**2,'l3^2':l3**2,'l2l3':l2*l3})
try_pslq("{1,pi2}", t2n, {'1':mpf(1),'pi2':pi**2})
try_pslq("{1,pi2,ln2,ln3}", t2n, {'1':mpf(1),'pi2':pi**2,'ln2':l2,'ln3':l3})
try_pslq("{1,pi2,ln2,ln3,ln2^2,ln3^2,ln2ln3}", t2n, {'1':mpf(1),'pi2':pi**2,'ln2':l2,'ln3':l3,'l2^2':l2**2,'l3^2':l3**2,'l2l3':l2*l3})
