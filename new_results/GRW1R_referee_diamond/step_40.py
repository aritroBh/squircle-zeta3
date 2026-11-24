import mpmath as mp, json
mp.mp.dps = 80
from fractions import Fraction
from math import factorial, comb
N=180
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]
r=[f450[n]/mp.factorial(n)-Hmp[n] for n in range(N+1)]
beta=[mp.mpf(0)]*(N+1); beta[0]=r[0]
twopi=2*mp.pi
for n in range(1,N+1):
    beta[n]=sum(r[n-2*j]*(twopi**(n-2*j))*comb(n-j-1,j) for j in range(0,(n-1)//2+1) if n-2*j>=1)
def w_of_z(z): return z/(mp.pi+mp.sqrt(mp.pi**2+z**2))
def R_map(z):
    w=w_of_z(z); tot=mp.mpf(0); pw=mp.mpf(1)
    for n in range(N+1): tot+=beta[n]*pw; pw*=w
    return tot
def Hz(z): return mp.sqrt(2)*mp.sqrt(1+mp.e**(-2*z))
def B_map(z): return Hz(z)+R_map(z)
def Hj_def(x, j):
    tot=mp.mpf(0)
    for m in range(0,N+1): tot+=Hmp[m]*x**(m+j)/mp.factorial(m+j)
    return tot
def ck(k):
    if k==0: return mp.mpf(1)
    rr=mp.mpf(1)
    for i in range(k): rr*=mp.mpf('0.5')-i
    return rr/mp.factorial(k)
CKS=[ck(k) for k in range(800)]
def Hj_neg_bessel(x, j, Kmax=700):
    ax=abs(x)
    tot=ax**(-mp.mpf(j)/2)*mp.besseli(j,2*mp.sqrt(ax))
    ssum=mp.mpf(0)
    for k in range(1,Kmax+1):
        y=(2*k-1)*ax
        ssum+=CKS[k]*y**(-mp.mpf(j)/2)*mp.besselj(j,2*mp.sqrt(y))
    return mp.sqrt(2)*x**j*(tot+ssum)
print("=== Lemma 6 identity check: defining series vs Bessel-I form, x<0 ===")
for x in [-1, -5, -20, -60, -85]:
    for j in [0,1]:
        v1=Hj_def(mp.mpf(x),j); v2=Hj_neg_bessel(mp.mpf(x),j)
        print(f"x={x}, j={j}: def={mp.nstr(v1,12)} bessel={mp.nstr(v2,12)} diff={mp.nstr(abs(v1-v2),3)}")
def kjfun(j,u):
    L=-mp.log(1+mp.e**(-u))
    return (1+mp.e**(-u))**(-2)*L**j/mp.factorial(j)
def BB(j, z, U):
    f=lambda u: mp.e**(-u)*kjfun(j,u)*(u**(-j) if j>0 else 1)*Hj_def(z*u, j)
    return mp.quad(f, [0, 1, 4, 10, U])
print("\n=== Lemma 7: B = sum BB_j, disk check ===")
for z in [mp.mpf('1.0'), mp.mpf('-1.0')]:
    direct=sum(f450[n]*z**n/mp.factorial(n) for n in range(N+1))
    approx=sum(BB(j,z,mp.mpf(20)) for j in range(0,8))
    print(f"z={z}: direct={mp.nstr(direct,15)} sumBB={mp.nstr(approx,15)} diff={mp.nstr(abs(direct-approx),3)}")
print("\n=== off-disk z=+-4 vs conformal map ===")
for z in [mp.mpf('4'), mp.mpf('-4')]:
    approx=sum(BB(j,z,mp.mpf(22)) for j in range(0,14))
    bm=B_map(z)
    print(f"z={z}: map={mp.nstr(bm,15)} sumBB={mp.nstr(approx,15)} diff={mp.nstr(abs(bm-approx),3)}")
with open('progress.log','a') as f:
    f.write("step_19: Lemma 6 Bessel-I form + Lemma 7 BB_j decomposition validated (disk & off-disk vs map)\n")
