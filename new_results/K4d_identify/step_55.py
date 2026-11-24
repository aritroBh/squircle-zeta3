import mpmath as mp, json
from mpmath import mpf
from fractions import Fraction as Fr
mp.mp.dps=75
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
def M(fr): return mpf(fr.numerator)/fr.denominator
z3=mp.zeta(3); pi=mp.pi; l2=mp.log(2); l3=mp.log(3); sq5=mp.sqrt(5); l43=mp.log(mpf(4)/3)
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
T192=getT('h192'); t2=T192[2]/sq5; t3=T192[3]/sq5
Li2m13=mp.polylog(2,mpf(-1)/3); Li3m13=mp.polylog(3,mpf(-1)/3)
def Aser(x):
    x=mpf(x); s=mpf(0); H=mpf(0); n=1
    while True:
        H+=mpf(1)/n; t=H*x**n/n**2; s+=t
        if abs(t)<mpf(10)**(-70) and n>8: break
        n+=1
    return s
Am13=Aser(mpf(-1)/3)
# VERIFY T2
t2_form = M(Fr(-254,675)) + M(Fr(79,8100))*pi**2 + M(Fr(-118,75))*Li2m13 + M(Fr(-4,15))*l43 + M(Fr(76,75))*l43**2
print("VERIFY T2/√5: |h192 - closed| =", mp.nstr(abs(t2-t2_form),3), " (closed: -254/675 +79pi2/8100 -118Li2(-1/3)/75 -4ln(4/3)/15 +76ln^2(4/3)/75)")
# T3 PSLQ: weight-3 basis, logs collapsed to ln(4/3); include Li2m13*l43 cross, pi2*l43, ln^3 43, z3, Li3m13, Am13
basis={'1':mpf(1),'z3':z3,'Li3m13':Li3m13,'Am13':Am13,'pi2l43':pi**2*l43,'Li2l43':Li2m13*l43,'l43^3':l43**3,'pi2':pi**2,'l43^2':l43**2,'l43':l43,'Li2m13':Li2m13}
names=list(basis)
def trypslq(tag,target,basis,tol,maxs=2000):
    nm=list(basis); x=[target]+[basis[n] for n in nm]
    rel,resid=pslq(x,tol=tol,maxsteps=maxs)
    if rel is None: print("   %s: none"%tag); return None
    h=max(abs(r) for r in rel); print("   %s: h=%.2e resid=%.1e"%(tag,float(h),float(resid)))
    return rel
print("\nT3/√5 =",mp.nstr(t3,50))
rel=trypslq("{w3,ln(4/3)} t3", t3, basis, mpf('1e-50'))
if rel:
    s=sum(rel[i+1]*basis[names[i]] for i in range(len(names)))
    print("   reconstruct diff:",mp.nstr(abs(-s/rel[0]-t3),3))
    for i,n in enumerate(names): print("     %-8s %s"%(n,Fr(rel[i+1],rel[0])))
