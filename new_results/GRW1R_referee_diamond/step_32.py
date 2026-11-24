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
rho=mp.pi/2; rhot=mp.mpf('2.9')*rho
def law2(n):
    return -mp.sqrt(2)*rho**(-n)*mp.mpf(n)**mp.mpf('-1.5')*mp.cos((2*n+1)*mp.pi/4)*(1 + (3-3*mp.pi*((-1)**n))/(8*n))
worst = mp.mpf(0); argworst=None
viol = []
for n in range(10, N+1):
    R = abs(Hmp[n]-law2(n))
    bound = mp.mpf('4.45')*rho**(-n)*mp.mpf(n)**mp.mpf('-3.5') + 135*rhot**(-n)
    ratio = R/bound
    if ratio > worst: worst, argworst = ratio, n
    if ratio > 1: viol.append(n)
print("claimed bound 4.45*u + 135*(2.9rho)^{-n}: worst ratio |R_n|/bound =", mp.nstr(worst,10), "at n =", argworst)
print("violations:", viol if viol else "NONE on 10..240")
# also check 7.06 bound without circle term
worst2 = mp.mpf(0); argw2=None
for n in range(10, N+1):
    R = abs(Hmp[n]-law2(n))
    ratio = R/(mp.mpf('7.06')*rho**(-n)*mp.mpf(n)**mp.mpf('-3.5'))
    if ratio > worst2: worst2, argw2 = ratio, n
print("7.06*u (no circle needed): worst ratio =", mp.nstr(worst2,10), "at n =", argw2)
with open('progress.log','a') as f:
    f.write(f"step_11: claimed 4.45+circle bound TRUE on 10..240 (worst {worst} at n={argworst}); 7.06 alone worst {worst2}\n")
