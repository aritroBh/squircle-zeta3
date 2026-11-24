import mpmath as mp, json
mp.mp.dps = 120
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
def etaD(n):
    if n==0: return mp.mpf('0.5')
    if n==1: return mp.log(2)
    return (1-mp.mpf(2)**(1-n))*mp.zeta(n)
def mid_series(z):
    return sum(Hmp[k]*(etaD(k)-1)*z**k for k in range(0, N+1))
mp.mp.dps = 60
def Hj_def(x,j):
    tot=mp.mpf(0)
    for m in range(0,N+1): tot+=Hmp[m]*x**(m+j)/mp.factorial(m+j)
    return tot
def kjfun(j,u):
    if u==0: return mp.mpf('0.25')*(-mp.log(2))**j/mp.factorial(j)
    L=-mp.log(1+mp.e**(-u))
    return (1+mp.e**(-u))**(-2)*L**j/mp.factorial(j)
def BB(j,z,U):
    f=lambda u: mp.e**(-u)*kjfun(j,u)*(u**(-j) if j>0 else 1)*Hj_def(z*u,j)
    return mp.quad(f, [0, 1, 4, 10, 18, U])
print("=== T2b: R_2 = sum_{j>=1} BB_j  vs  B_map - H - sum_k H_k(eta_D(k)-1) z^k ===")
for z in [mp.mpf('2.0'), mp.mpf('2.5'), mp.mpf('3.0')]:
    U = mp.mpf(28)
    R2q = sum(BB(j, z, U) for j in range(1, 16))
    rhs = B_map(z) - Hz(z) - mid_series(z)
    print(f"z={z}: R2(quad)={mp.nstr(R2q,14)}  map-side={mp.nstr(rhs,14)}  diff={mp.nstr(abs(R2q-rhs),3)}")
print("\n=== claim A spot check, eta=0.37 ===")
eta = mp.mpf('0.37')
def Kfun(e_,u): return (1+mp.e**(-u))**(-2-e_)
lap = mp.quad(lambda u: mp.e**(-u)*Hz(eta*u)*Kfun(eta,u), [0, 5, 15, mp.inf])
f = lambda t: mp.sqrt(t**(2*eta)+(1-t)**(2*eta))
P = 4*(1+eta)*mp.quad(f, [0, mp.mpf('0.5'), 1])
Iv = P/(4*mp.sqrt(2)*(1+eta))
print(f"Laplace: {mp.nstr(lap,25)}\nperim:   {mp.nstr(Iv,25)}\ndiff: {mp.nstr(abs(lap-Iv),3)}")
with open('progress.log','a') as f:
    f.write("step_28: T2b decomposition check at z=2,2.5,3 + claim A spot check done\n")
