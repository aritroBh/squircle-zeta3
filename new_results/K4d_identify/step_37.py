import mpmath as mp, json
from mpmath import mpf
mp.mp.dps = 60
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0, WORK)
from pslq_lib import pslq

def Afunc(x):
    # sum_{n>=1} H_n x^n / n^2
    x=mpf(x); s=mpf(0); H=mpf(0); n=1
    while True:
        H+=mpf(1)/n
        term=H*x**n/n**2
        s+=term
        if abs(term)<mpf(10)**(-mp.mp.dps+5) and n>5: break
        n+=1
    return s

# verify A(-1) = -5 zeta3/8 and A(1/2)= zeta3 - pi^2 ln2/12
z3=mp.zeta(3); pi=mp.pi; l2=mp.log(2); l3=mp.log(3)
print("check A(-1) =", mp.nstr(Afunc(-1),30), " vs -5z3/8=", mp.nstr(-5*z3/8,30))
print("check A(1/2)=", mp.nstr(Afunc(0.5),30), " vs z3-pi^2 l2/12=", mp.nstr(z3-pi**2*l2/12,30))

# T2/S0 from h96
def getT(fn):
    SA=[mpf(s) for s in json.load(open(WORK+'/'+fn))['Sdk']]
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
    u=[mpf(0),1,-1,1,-1]; Lu=[L*x for x in u]; e=[mpf(0)]*5; e[0]=1; pw=[mpf(0)]*5; pw[0]=1
    for k in range(1,5):
        pw=sermul(pw,Lu,4)
        for i in range(4,-1,-1): e[i]+=pw[i]/mp.factorial(k)
    return sermul(list(SA),serinv(e,4),4)
T96=getT('modelA_h96.json'); S0=T96[0]
t2=T96[2]/S0
Li2_m13=mp.polylog(2,mpf(-1)/3)
consts={'1':mpf(1),'pi2':pi**2,'Li2(-1/3)':Li2_m13,'ln2^2':l2**2,'ln3^2':l3**2,'ln2ln3':l2*l3}
print("\nt2 =", mp.nstr(t2,45))
for nm,c in consts.items(): print("  %s = %s"%(nm,mp.nstr(c,30)))
# PSLQ [t2, consts]
x=[t2]+[consts[k] for k in consts]
rel,resid=pslq(x, tol=mpf('1e-38'), maxsteps=2000)
print("\nPSLQ(t2, 1, pi2, Li2(-1/3), ln2^2, ln3^2, ln2ln3):")
print("  relation:", rel)
print("  residual:", mp.nstr(resid,3) if resid else None)
