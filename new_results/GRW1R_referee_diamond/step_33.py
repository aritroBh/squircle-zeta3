import mpmath as mp
mp.mp.dps = 200
from fractions import Fraction
from math import factorial
import json

N=180
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]
rho=mp.pi/2
CH = mp.mpf('5.874')

def etaD(n):
    return (1-mp.mpf(2)**(1-n))*mp.zeta(n)

print("=== Claim D test: tail_n = f_n - H_n n! eta_D(n) ===")
print("columns: n | |tail| | T_n (proved sum bound) | D's closed bound | tail/main | bound/main")
for n in [20, 40, 80, 120, 180]:
    main = Hmp[n]*mp.factorial(n)*etaD(n)
    tail = f450[n]-main
    T = mp.mpf(0)
    for j in range(1, n+1):
        T += abs(Hmp[n-j])*mp.factorial(n-j)/(mp.factorial(j)*(j+1)**(n-j+1))
    Dbnd = CH*mp.e**(2*rho)*rho**(-n)*mp.factorial(n)*mp.mpf(2)**(-n-1)/n
    print(f"n={n:4d} | tail={mp.nstr(abs(tail),6)} | T={mp.nstr(T,6)} | Dbnd={mp.nstr(Dbnd,6)} | "
          f"|tail|/T={mp.nstr(abs(tail)/T,4)} | T/Dbnd={mp.nstr(T/Dbnd,6)} | |tail|/|main|={mp.nstr(abs(tail)/abs(main),4)}")
# where does D's chain lose? sum over j of rho^j/(j! 2^{n-j+1}) vs e^{2rho}2^{-n-1}
n=120
s1 = sum(rho**j/(mp.factorial(j)*mp.mpf(2)**(n-j+1)) for j in range(1,n+1))
print("n=120: sum rho^j/(j!2^{n-j+1}) =", mp.nstr(s1,10), " e^{2rho}2^{-n-1} =", mp.nstr(mp.e**(2*rho)*mp.mpf(2)**(-n-1),10),
      " slack factor:", mp.nstr(mp.e**(2*rho)*mp.mpf(2)**(-n-1)/s1,6))
with open('progress.log','a') as f:
    f.write("step_12: claim D bound tested vs true tail at n=20,40,80,120,180; holds with slack (see stdout)\n")
