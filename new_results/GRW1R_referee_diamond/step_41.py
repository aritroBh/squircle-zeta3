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
def ck(k):
    if k==0: return mp.mpf(1)
    rr=mp.mpf(1)
    for i in range(k): rr*=mp.mpf('0.5')-i
    return rr/mp.factorial(k)
CKS=[ck(k) for k in range(900)]
def Hj_def(x,j):
    tot=mp.mpf(0)
    for m in range(0,N+1): tot+=Hmp[m]*x**(m+j)/mp.factorial(m+j)
    return tot
def Hj_neg_bessel(x,j,Kmax=800):
    ax=abs(x)
    tot=ax**(-mp.mpf(j)/2)*mp.besseli(j,2*mp.sqrt(ax))
    ssum=mp.mpf(0)
    for k in range(1,Kmax+1):
        y=(2*k-1)*ax
        ssum+=CKS[k]*y**(-mp.mpf(j)/2)*mp.besselj(j,2*mp.sqrt(y))
    return mp.sqrt(2)*x**j*(tot+ssum)
def Hj(x,j):
    if x >= -85: return Hj_def(x,j)
    return Hj_neg_bessel(x,j)
def kjfun(j,u):
    L=-mp.log(1+mp.e**(-u))
    return (1+mp.e**(-u))**(-2)*L**j/mp.factorial(j)
def BB(j,z,U):
    f=lambda u: mp.e**(-u)*kjfun(j,u)*(u**(-j) if j>0 else 1)*Hj(z*u,j)
    return mp.quad(f, [0, 1, 4, 10, 20, U])
print("=== redo z=-4 with U=44 ===")
z=mp.mpf('-4')
approx=sum(BB(j,z,mp.mpf(44)) for j in range(0,14))
bm=B_map(z)
print(f"z=-4: map={mp.nstr(bm,16)}  sumBB={mp.nstr(approx,16)}  diff={mp.nstr(abs(bm-approx),3)}")
print("=== redo z=-1 disk check with U=32 ===")
z=mp.mpf('-1')
direct=sum(f450[n]*z**n/mp.factorial(n) for n in range(N+1))
approx=sum(BB(j,z,mp.mpf(32)) for j in range(0,8))
print(f"z=-1: direct={mp.nstr(direct,15)} sumBB={mp.nstr(approx,15)} diff={mp.nstr(abs(direct-approx),3)}")
with open('progress.log','a') as f:
    f.write("step_20: off-disk z=-4 BB-sum vs map fixed (U=44, split Hj); disk z=-1 redone U=32\n")
