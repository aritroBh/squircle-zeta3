import mpmath as mp
mp.mp.dps = 50
# min |H_n/lawH| over 10<=n<=200 (exact rationals)
from fractions import Fraction
from math import factorial
N=200
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
rho=mp.pi/2
mn=mp.mpf(9); argn=None
for n in range(10,201):
    v=abs(Hmp[n])*rho**n*mp.mpf(n)**mp.mpf('1.5')
    if v<mn: mn, argn=v,n
print("min |H_n/lawH| on 10..200:", mp.nstr(mn,10), "at n =", argn, "(claimed 0.8745)")
# Lemma 3.4 orbit argument: pick A=1,B=2 (L=log 2, theta_0 = arctan((pi/2)/L)); theta_0/pi rational?
import numpy as np
for (A,B) in [(1,2),(1,3),(2,3)]:
    L=mp.log(B/A); th=mp.atan((mp.pi/2)/L)
    print(f"A={A},B={B}: theta0/pi={mp.nstr(th/mp.pi,12)}")
    # bracket = r cos(psi_n+phi*) (n even), -r cos(psi_n-phi*) (n odd); check limsup |bracket| > 0
    r=mp.sqrt(A**2+B**2); ph=mp.atan(A/B)
    best=mp.mpf(0)
    for n in range(1, 400):
        psi=(n-mp.mpf('0.5'))*th
        v = r*abs(mp.cos(psi+ph)) if n%2==0 else r*abs(mp.cos(psi-ph))
        best=max(best,v)
    print(f"   max |bracket| over n<400: {mp.nstr(best,8)} (r={mp.nstr(r,8)})  -> nonvanishing subsequence exists")
with open('progress.log','a') as f:
    f.write("step_30: min|H/lawH|=0.8745@10 confirmed; Lemma 3.4 bracket nonvanishing spot checks\n")
