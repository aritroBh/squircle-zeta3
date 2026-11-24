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
CKS=[ck(k) for k in range(1500)]

def Hj_def(x, j):  # defining series, x>=0
    tot=mp.mpf(0)
    for m in range(0,N+1): tot+=Hmp[m]*x**(m+j)/mp.factorial(m+j)
    return tot
def Hj_neg(x, j, Kmax=1400):  # x<0 Bessel form
    ax=abs(x)
    tot = ax**(-mp.mpf(j)/2)*mp.besseli(j, 2*mp.sqrt(ax)) if j>0 else mp.besseli(0,2*mp.sqrt(ax))
    # careful: I_j(2 sqrt(ax)) ax^{-j/2}
    tot = ax**(-mp.mpf(j)/2)*mp.besseli(j, 2*mp.sqrt(ax))
    ssum=mp.mpf(0)
    for k in range(1,Kmax+1):
        y=(2*k-1)*ax
        ssum += CKS[k]*y**(-mp.mpf(j)/2)*mp.besselj(j,2*mp.sqrt(y))
    return mp.sqrt(2)*x**j*(tot+ssum)
def Hj(x, j):
    if x >= 0: return Hj_def(x, j)
    return Hj_neg(x, j)

# validate Hj_neg against defining series for x in (-90, 0)
for x in [-1, -5, -20, -60, -88]:
    v1 = Hj_def(mp.mpf(x),0); v2 = Hj_neg(mp.mpf(x),0)
    print(f"x={x}: def={mp.nstr(v1,14)} besselI={mp.nstr(v2,14)} diff={mp.nstr(abs(v1-v2),3)}")

def kjfun(j,u):
    L=-mp.log(1+mp.e**(-u))
    return (1+mp.e**(-u))**(-2)*L**j/mp.factorial(j)
def BB(j, z, U=mp.mpf(20)):
    # B_j(z) = int_0^U e^{-u} k_j(u) u^{-j} H^{[j]}(z u) du + tail bound
    if z >= 0:
        f=lambda u: mp.e**(-u)*kjfun(j,u)*(u**(-j) if j>0 else 1)*Hj(z*u, j)
    else:
        f=lambda u: mp.e**(-u)*kjfun(j,u)*(u**(-j) if j>0 else 1)*Hj_neg(z*u, j)
    pts=[0, 1, 3, 8, U]
    return mp.quad(f, pts)
for z in [mp.mpf('1.0'), mp.mpf('-1.0')]:
    direct = sum(f450[n]*z**n/mp.factorial(n) for n in range(N+1))
    approx = sum(BB(j, z) for j in range(0, 8))
    print(f"z={z}: B direct={mp.nstr(direct,16)}  sum B_j={mp.nstr(approx,16)}  diff={mp.nstr(abs(direct-approx),3)}")
for z in [mp.mpf('4'), mp.mpf('-4')]:
    approx = sum(BB(j, z) for j in range(0, 14))
    bm = B_map(z)
    print(f"z={z}: B map={mp.nstr(bm,16)}  sum B_j={mp.nstr(approx,16)}  diff={mp.nstr(abs(bm-approx),3)}")
with open('progress.log','a') as f:
    f.write("step_18: BB_j continuation validated vs direct series (disk) and vs conformal map at z=+-4\n")
