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
Li2m13=mp.polylog(2,mpf(-1)/3)
basis={'1':mpf(1),'pi2':pi**2,'Li2m13':Li2m13,'ln2':l2,'ln3':l3,'l2^2':l2**2,'l3^2':l3**2,'l2l3':l2*l3}
names=list(basis)
def run(tag):
    T=getT(tag); t2=T[2]/sq5
    x=[t2]+[basis[n] for n in names]
    rel,resid=pslq(x,tol=mpf('1e-52'),maxsteps=1500)
    return t2,rel,resid
rels={}
for tag in ['h96','h128','h192']:
    t2,rel,resid=run(tag)
    rels[tag]=rel
    print("%s: resid=%.1e  rel=%s"%(tag,float(resid),rel))
# Check if h192 relation reproduces t2 exactly
T=getT('h192'); t2=T[2]/sq5
rel=rels['h192']
if rel:
    s=sum(rel[i+1]*basis[names[i]] for i in range(len(names)))
    print("\nreconstruct vs t2: diff=",mp.nstr(abs(-s/rel[0]-t2),3))
    # simplify: rationalize small coeffs
    print("coeffs rel[0]=%d"%rel[0])
    for i,n in enumerate(names):
        from fractions import Fraction
        print("   %-7s %s"%(n, Fraction(rel[i+1],rel[0])))
