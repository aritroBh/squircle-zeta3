import mpmath as mp
mp.mp.dps = 60
def ck(k):
    if k==0: return mp.mpf(1)
    r = mp.mpf(1)
    for i in range(k): r *= mp.mpf('0.5')-i
    return r/mp.factorial(k)
cks = [ck(k) for k in range(6000)]

def Hbessel_pos_rich(x, j=0, Ks=(2000,3000,4000,5000)):
    # partial sums at several K, Richardson-extrapolate K->inf assuming error ~ K^{-3/4}
    vals = []
    for K in Ks:
        tot = mp.mpf(1)/mp.factorial(j)
        for k in range(1, K+1):
            y = 2*k*x
            tot += cks[k]*y**(-mp.mpf(j)/2)*mp.besselj(j, 2*mp.sqrt(y))
        vals.append(mp.sqrt(2)*x**j*tot)
    # fit error = a K^{-3/4} + b K^{-1}: use last 3
    import numpy as np
    # simple: extrapolate with model v(K) = V + a K^{-3/4} + b K^{-1} using Ks[-3:]
    M = mp.matrix(3,3); b = mp.matrix(3,1)
    for i,K in enumerate(Ks[-3:]):
        M[i,0]=1; M[i,1]=mp.mpf(K)**mp.mpf('-0.75'); M[i,2]=mp.mpf(K)**-1
        b[i]=vals[-3:][i]
    sol = mp.lu_solve(M,b)
    return sol[0], vals

# validate at x=80 against defining series
mp.mp.dps = 100
from fractions import Fraction
from math import factorial
N=180
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
def Hpow(x, j=0):
    tot=mp.mpf(0)
    for m in range(0,N+1): tot += Hmp[m]*x**(m+j)/mp.factorial(m+j)
    return tot
mp.mp.dps = 60
for x in [40, 80]:
    v1 = Hpow(mp.mpf(x),0)
    v2, vals = Hbessel_pos_rich(mp.mpf(x),0)
    print(f"x={x}: defining={mp.nstr(v1,20)}")
    print(f"     rich={mp.nstr(v2,20)}  diff={mp.nstr(abs(v1-v2),4)}  partials={[mp.nstr(v,8) for v in vals]}")
# large x sanity: should stay <= 2 sqrt2 = 2.828
for x in [200, 400, 800, 1200]:
    v2, _ = Hbessel_pos_rich(mp.mpf(x),0)
    print(f"x={x}: H^[0] ~ {mp.nstr(v2,10)}   (bound 2.8284)")
with open('progress.log','a') as f:
    f.write("step_15: Richardson-accelerated Bessel series for H^[0] on R+, large-x sanity vs 2sqrt2 bound\n")
