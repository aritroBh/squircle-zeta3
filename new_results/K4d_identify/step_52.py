import mpmath as mp, json
from mpmath import mpf
mp.mp.dps=75
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
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
def getT(tag):
    SA=[mpf(s) for s in json.load(open(WORK+'/modelA_%s.json'%tag))['Sdk']]
    return sermul(list(SA),finv,4)
T192=getT('h192'); T128=getT('h128')
t2=T192[2]/sq5; t3=T192[3]/sq5
print("t2/sqrt5(h192) =",mp.nstr(t2,60))
print("t3/sqrt5(h192) =",mp.nstr(t3,60))
def try_pslq(tag,target,basis,tol=mpf('1e-52')):
    names=list(basis); x=[target]+[basis[n] for n in names]
    rel,resid=pslq(x,tol=tol,maxsteps=1000)
    if rel is None: print("   %s: none"%tag); return None
    h=max(abs(r) for r in rel)
    mark=' <<<' if h<mpf('1e6') else ''
    print("   %s: h=%.2e resid=%.1e%s"%(tag,float(h),float(resid),mark))
    return rel,resid,h
# retest the borderline one
bb={'1':mpf(1),'pi2':pi**2,'ln2':l2,'ln3':l3,'l2^2':l2**2,'l3^2':l3**2,'l2l3':l2*l3}
r=try_pslq("{1,pi2,ln2,ln3,ln2^2,ln3^2,ln2ln3} t2", t2, bb)
# test with Li2(-1/3)
Li2m13=mp.polylog(2,mpf(-1)/3)
try_pslq("{1,pi2,Li2(-1/3),ln2,ln3,ln2^2,ln3^2,ln2ln3} t2", t2, dict(bb, **{'Li2m13':Li2m13}))
try_pslq("{1,pi2,Li2(-1/3)} t2", t2, {'1':mpf(1),'pi2':pi**2,'Li2m13':Li2m13})
try_pslq("{1,pi2,Li2(-1/3),ln2^2,ln3^2,ln2ln3} t2", t2, {'1':mpf(1),'pi2':pi**2,'Li2m13':Li2m13,'l2^2':l2**2,'l3^2':l3**2,'l2l3':l2*l3})
