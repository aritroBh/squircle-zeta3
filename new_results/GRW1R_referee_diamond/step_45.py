import mpmath as mp, json
from fractions import Fraction
from math import factorial
mp.mp.dps = 60
N=30
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]
def kjfun(j,u):
    if u==0: return mp.mpf('0.25')*(-mp.log(2))**j/mp.factorial(j)
    L=-mp.log(1+mp.e**(-u))
    return (1+mp.e**(-u))**(-2)*L**j/mp.factorial(j)
def M(n,j):
    f=lambda u: mp.e**(-u)*u**(n-j)*kjfun(j,u)
    v, err = mp.quad(f, [0, mp.inf], error=True, maxdegree=8)
    return v, err
worst = mp.mpf(0)
for n in [1,2,3,5,8]:
    tot=mp.mpf(0); terr=mp.mpf(0)
    for j in range(0,n+1):
        v,err = M(n,j)
        tot += Hmp[n-j]*v; terr += abs(Hmp[n-j])*err
    rel = abs(tot-f450[n])/abs(f450[n])
    worst = max(worst, rel)
    print(f"n={n}: rel diff={mp.nstr(rel,3)}  (quad err est {mp.nstr(terr/abs(f450[n]),3)})")
# exact spot values: M(1,1) = (ln2-1)/2, M(2,1)?
v,_=M(1,1); print("M(1,1) vs (ln2-1)/2:", mp.nstr(v,20), mp.nstr((mp.log(2)-1)/2,20), mp.nstr(abs(v-(mp.log(2)-1)/2),3))
with open('progress.log','a') as f:
    f.write(f"step_24: Lemma 2 own-quadrature check n=1..8 worst rel {mp.nstr(worst,3)}; M(1,1) exact match\n")
