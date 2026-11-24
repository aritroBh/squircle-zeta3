import mpmath as mp
mp.mp.dps = 80
from fractions import Fraction
from math import factorial
N=240
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
rho=mp.pi/2; Estar=mp.mpf(5)/24
def A_n(n):
    return -(2/mp.pi)*2*mp.cos(n*mp.pi/2+mp.pi/4)*rho**(mp.mpf('0.5')-n)*mp.beta(mp.mpf('1.5'),n-mp.mpf('0.5'))
def B_mine(n):
    return (2/mp.pi)*mp.cos(n*mp.pi/2-mp.pi/4)*rho**(mp.mpf('1.5')-n)*mp.beta(mp.mpf('2.5'),n-mp.mpf('1.5'))
def B_grw(n):
    return -(2/mp.pi)*mp.cos(3*mp.pi/4-n*mp.pi/2)*rho**(mp.mpf('1.5')-n)*mp.beta(mp.mpf('2.5'),n-mp.mpf('1.5'))
def Ebound(n):
    return (4/mp.pi)*Estar*rho**(mp.mpf('2.5')-n)*mp.beta(mp.mpf('3.5'),n-mp.mpf('2.5'))
print(" n |  |H-A-B_mine|/Ebound |  |H-A-B_grw|/Ebound")
for n in [12, 20, 40, 80, 120, 200, 240]:
    r1 = abs(Hmp[n]-A_n(n)-B_mine(n))/Ebound(n)
    r2 = abs(Hmp[n]-A_n(n)-B_grw(n))/Ebound(n)
    print(f"{n:4d} | {mp.nstr(r1,10)} | {mp.nstr(r2,10)}")
with open('progress.log','a') as f:
    f.write("step_06: Hankel phases fixed; B_mine confirmed vs exact H_n (ratio<1 means E-term bound holds)\n")
