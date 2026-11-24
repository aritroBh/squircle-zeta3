import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction
from math import factorial
N = 30
KK = N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
rho=mp.pi/2
def law1(n):
    return -mp.sqrt(2)*rho**(-n)*mp.mpf(n)**mp.mpf('-1.5')*mp.cos((2*n+1)*mp.pi/4)
print(" n   H_n                 law1                H_n/law1")
for n in range(1, 21):
    print(f"{n:3d}  {mp.nstr(Hmp[n],12):22s} {mp.nstr(law1(n),12):22s} {mp.nstr(Hmp[n]/law1(n),10)}")
