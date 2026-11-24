import mpmath as mp
mp.mp.dps = 80
from fractions import Fraction
from math import factorial
N = 240
KK = N+2
a = [Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
rho=mp.pi/2; rhot=mp.mpf('2.9')*rho; Estar=mp.mpf(5)/24

def Ainf(n):
    ph = -2*mp.cos(n*mp.pi/2 + mp.pi/4)      # 2Re[-e^{i pi/4} i^n]
    return (2/mp.pi)*ph*rho**(mp.mpf('0.5')-n)*mp.beta(mp.mpf('1.5'),n-mp.mpf('0.5'))
# fit B phase: B_n = (2/pi)*rho^{3/2-n} B(5/2,n-3/2) * phB(n); try phB = 2Re[lam i^n], lam in a few choices
cands = {
 'i*e^{ipi/4}': lambda n: 2*mp.cos(n*mp.pi/2+mp.pi/2+mp.pi/4),
 '-i*e^{ipi/4}': lambda n: -2*mp.cos(n*mp.pi/2+mp.pi/2+mp.pi/4),
 'e^{ipi/4}': lambda n: 2*mp.cos(n*mp.pi/2+mp.pi/4),
 '-e^{ipi/4}': lambda n: -2*mp.cos(n*mp.pi/2+mp.pi/4),
}
def Eterm_true_bound(n):
    return (4/mp.pi)*Estar*rho**(mp.mpf('2.5')-n)*mp.beta(mp.mpf('3.5'),n-mp.mpf('2.5'))
# determine phB by requiring |H_n - A - B| <= Ebound + tiny tails for several n
for name, ph in cands.items():
    ok = True
    worst = mp.mpf(0)
    for n in [12, 20, 40, 80]:
        B = (2/mp.pi)*ph(n)*rho**(mp.mpf('1.5')-n)*mp.beta(mp.mpf('2.5'),n-mp.mpf('1.5'))
        resid = abs(Hmp[n]-Ainf(n)-B)
        ratio = resid/Eterm_true_bound(n)
        worst = max(worst, ratio)
    print(name, " worst |H-A-B|/Ebound =", mp.nstr(worst,6))
