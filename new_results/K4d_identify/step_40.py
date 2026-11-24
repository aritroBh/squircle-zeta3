import mpmath as mp, json, time
from mpmath import mpf
mp.mp.dps=65
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
z3=mp.zeta(3); pi=mp.pi; l2=mp.log(2); l3=mp.log(3)
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
T=sermul(list(SA),finv,4); S0=T[0]
t2=T[2]/S0; t3=T[3]/S0
print("t2 =",mp.nstr(t2,55))
print("t3 =",mp.nstr(t3,55))
Li2m13=mp.polylog(2,mpf(-1)/3)
basis={'1':mpf(1),'pi2':pi**2,'Li2(-1/3)':Li2m13,'ln2^2':l2**2,'ln3^2':l3**2,'ln2ln3':l2*l3}
names=list(basis)
x=[t2]+[basis[n] for n in names]
t0=time.time()
rel,resid=pslq(x,tol=mpf('1e-46'),maxsteps=500)
print("\nPSLQ t2 on [1,pi2,Li2(-1/3),ln2^2,ln3^2,ln2ln3] (%.1fs):"%(time.time()-t0))
print("  rel:",rel)
print("  resid:",mp.nstr(resid,3) if resid else None)
if rel:
    # reconstruct t2
    s=sum(rel[i+1]*basis[names[i]] for i in range(len(names)))
    print("  check: -(%s)/%d = %s"%( ' + '.join('%d*%s'%(rel[i+1],names[i]) for i in range(len(names))), rel[0], mp.nstr(-s/rel[0],50)))
    print("  t2    =",mp.nstr(t2,50))
    print("  diff  =",mp.nstr(abs(-s/rel[0]-t2),3))
