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
rho=mp.pi/2; rhot=mp.mpf('2.9')*rho; Estar=mp.mpf(5)/24
def A_n(n):
    return -(2/mp.pi)*2*mp.cos(n*mp.pi/2+mp.pi/4)*rho**(mp.mpf('0.5')-n)*mp.beta(mp.mpf('1.5'),n-mp.mpf('0.5'))
def A_grw(n):  # (2/pi) 2Re[e^{i pi/4} i^{-n}] rho^{1/2-n} B(3/2,n-1/2)
    return (2/mp.pi)*2*mp.cos(mp.pi/4-n*mp.pi/2)*rho**(mp.mpf('0.5')-n)*mp.beta(mp.mpf('1.5'),n-mp.mpf('0.5'))
def B_mine(n):
    return (2/mp.pi)*mp.cos(n*mp.pi/2-mp.pi/4)*rho**(mp.mpf('1.5')-n)*mp.beta(mp.mpf('2.5'),n-mp.mpf('1.5'))
def Ebound(n):
    return (4/mp.pi)*Estar*rho**(mp.mpf('2.5')-n)*mp.beta(mp.mpf('3.5'),n-mp.mpf('2.5'))
# check A_grw at odd n
for n in [11, 13]:
    print(n, "A_grw vs H: resid/Ebound =", mp.nstr(abs(Hmp[n]-A_grw(n)-B_mine(n))/Ebound(n),8),
          "  A_mine:", mp.nstr(abs(Hmp[n]-A_n(n)-B_mine(n))/Ebound(n),8))
print("A_grw(11), A_mine(11), H_11:", mp.nstr(A_grw(11),12), mp.nstr(A_n(11),12), mp.nstr(Hmp[11],12))
# now honest remainder
def law2(n):
    return -mp.sqrt(2)*rho**(-n)*mp.mpf(n)**mp.mpf('-1.5')*mp.cos((2*n+1)*mp.pi/4)*(1 + (3-3*mp.pi*((-1)**n))/(8*n))
def tailA(n): return (4/mp.pi)*rhot**(mp.mpf('0.5')-n)/(n-mp.mpf('0.5'))
def tailB(n): return (2/mp.pi)*rhot**(mp.mpf('1.5')-n)/(n-mp.mpf('1.5'))
def circ(n): return 135*rhot**(-n)
print("\n n | |A+B-law2|/u | Eb/u | tails/u | honest/u | true|R|/u | u=rho^-n n^-3.5")
for n in [10,11,12,13,15,20,30,50,80,120,180,240]:
    u = rho**(-n)*mp.mpf(n)**mp.mpf('-3.5')
    ABm = abs(A_n(n)+B_mine(n)-law2(n))/u
    Eb = Ebound(n)/u
    tl = (tailA(n)+tailB(n))/u
    honest = ABm+Eb+tl
    trueR = abs(Hmp[n]-law2(n))/u
    print(f"{n:4d} | {mp.nstr(ABm,9)} | {mp.nstr(Eb,9)} | {mp.nstr(tl,5)} | {mp.nstr(honest,9)} | {mp.nstr(trueR,9)}")
with open('progress.log','a') as f:
    f.write("step_07: honest remainder assembled: |A+B-law2|+Ebound+tails; see stdout\n")
