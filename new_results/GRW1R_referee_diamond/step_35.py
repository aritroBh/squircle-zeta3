import mpmath as mp
mp.mp.dps = 80
from fractions import Fraction
from math import factorial
N=180
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]

def Hpow(x, j=0, N=N):  # defining series H^{[j]}(x) = sum H_m x^{m+j}/(m+j)!
    tot = mp.mpf(0)
    for m in range(N, -1, -1):
        tot = tot*x + Hmp[m]/mp.factorial(m+j) if False else tot
    # direct sum
    tot = mp.mpf(0)
    for m in range(0, N+1):
        tot += Hmp[m]*x**(m+j)/mp.factorial(m+j)
    return tot

def ck(k):  # binom(1/2,k)
    if k==0: return mp.mpf(1)
    r = mp.mpf(1)
    for i in range(k):
        r *= mp.mpf('0.5')-i
    return r/mp.factorial(k)

def Hbessel_pos(x, j=0, Kmax=4000):
    # H^{[j]}(x) = sqrt2 x^j sum_{k>=0} c_k F_j(-2 k x), F_j(y<0)=|y|^{-j/2} J_j(2 sqrt|y|); k=0 term = 1/j!
    tot = mp.mpf(1)/mp.factorial(j)
    for k in range(1, Kmax+1):
        y = 2*k*x
        F = y**(-mp.mpf(j)/2)*mp.besselj(j, 2*mp.sqrt(y))
        tot += ck(k)*F
    return mp.sqrt(2)*x**j*tot

for x in [5, 20, 40, 80]:
    v1 = Hpow(mp.mpf(x),0)
    v2 = Hbessel_pos(mp.mpf(x),0)
    print(f"x={x}: defining={mp.nstr(v1,18)}  bessel={mp.nstr(v2,18)}  diff={mp.nstr(abs(v1-v2),4)}")
with open('progress.log','a') as f:
    f.write("step_14: H^{[0]} defining series vs Bessel series on R+; see stdout\n")
