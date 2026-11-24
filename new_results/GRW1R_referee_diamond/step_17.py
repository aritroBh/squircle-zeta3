import mpmath as mp
mp.mp.dps = 80
from fractions import Fraction
from math import factorial

# exact H_n (rationals)
N = 220
KK = N+2
a = [Fraction(0)]*KK; a[0] = Fraction(1)
for k in range(1, KK): a[k] = Fraction((-2)**k, 2*factorial(k))
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK): s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hmp = [mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]

rho = mp.pi/2
L = mp.mpf('1.9')*rho
rhot = mp.mpf('2.9')*rho

def law2(n):
    c = mp.cos((2*n+1)*mp.pi/4)
    return -mp.sqrt(2)*rho**(-n)*mp.mpf(n)**mp.mpf('-1.5')*c*(1 + (3-3*mp.pi*((-1)**n))/(8*n))

# Hankel main pieces with full Beta integrals:
# H_n = A_n^inf + B_n^inf + (E-term) + (tail corrections) + circle
# A_n^inf = (2/pi)*2Re[e^{i pi/4} i^{-n}] * rho^{1/2-n} * B(3/2, n-1/2)   [GRW1 convention i^{-n}]
# B_n^inf = (2/pi)*2Re[-(i/2) e^{i pi/4} i^{-n}] ... hmm GRW1: -(1/pi) rho^{3/2-n} B(5/2,n-3/2) 2Re[i e^{i pi/4} i^{-n}]
def Ainf(n):
    ph = 2*mp.cos(mp.pi/4 + n*mp.pi/2)   # 2Re[e^{i pi/4} i^{-n}]
    return (2/mp.pi)*ph*rho**(mp.mpf('0.5')-n)*mp.beta(mp.mpf('1.5'), n-mp.mpf('0.5'))
def Binf(n):
    ph = 2*mp.cos(mp.pi/2 + mp.pi/4 + n*mp.pi/2)  # 2Re[i e^{i pi/4} i^{-n}]
    return -(1/mp.pi)*ph*rho**(mp.mpf('1.5')-n)*mp.beta(mp.mpf('2.5'), n-mp.mpf('1.5'))

print("check A+B vs H_n - (E stuff): n=12")
n=12
print(mp.nstr(Hmp[n],30), mp.nstr(Ainf(n)+Binf(n),30))
print("true R_n = H_n - law2:", mp.nstr(Hmp[n]-law2(n),20))
print("A+B-law2:", mp.nstr(Ainf(n)+Binf(n)-law2(n),20))

Estar = mp.mpf(5)/24
def Ebound(n):
    return (4/mp.pi)*Estar*rho**(mp.mpf('2.5')-n)*mp.beta(mp.mpf('3.5'), n-mp.mpf('2.5'))
def tailA(n):  # (4/pi) int_L^inf sqrt(s) (rho+s)^{-n-1} ds <= (4/pi)(2.9rho)^{1/2-n}/(n-1/2)
    return (4/mp.pi)*rhot**(mp.mpf('0.5')-n)/(n-mp.mpf('0.5'))
def tailB(n):  # (2/pi) int_L^inf s^{3/2}(rho+s)^{-n-1} ds <= (2/pi)(2.9rho)^{3/2-n}/(n-3/2)
    return (2/mp.pi)*rhot**(mp.mpf('1.5')-n)/(n-mp.mpf('1.5'))
def circ(n):
    return 135*rhot**(-n)

print("\nn | true|R|/(rho^-n n^-3.5) | Ebound | A+B-law2 piece | tails | honest total | claimed 4.45 |")
for n in [10, 11, 12, 15, 20, 30, 50, 80, 120, 180, 220]:
    unit = rho**(-n)*mp.mpf(n)**mp.mpf('-3.5')
    trueR = abs(Hmp[n]-law2(n))/unit
    Eb = Ebound(n)/unit
    ABm = abs(Ainf(n)+Binf(n)-law2(n))/unit
    tl = (tailA(n)+tailB(n))/unit
    honest = Eb + ABm + tl
    print(f"{n:4d} | {mp.nstr(trueR,10)} | {mp.nstr(Eb,10)} | {mp.nstr(ABm,10)} | {mp.nstr(tl,6)} | {mp.nstr(honest,10)} | 4.45")
with open('progress.log','a') as f:
    f.write("step_05: honest Darboux remainder decomposition computed; see stdout\n")
