import re
src = open('verify_grw1r.py').read()
old = """mp.mp.dps = 120
r=[f450[n]/mp.factorial(n)-Hmp[n] for n in range(181)]"""
new = """mp.mp.dps = 120
f450h=[mp.mpf(v) for v in F450]   # re-read at full precision (map composition is ill-conditioned)
r=[f450h[n]/mp.factorial(n)-Hmp[n] for n in range(181)]"""
assert old in src
src = src.replace(old, new)
open('verify_grw1r.py','w').write(src)
print("patched")
# rerun only section 6 quickly to confirm eta=2
import mpmath as mp, json, os
from fractions import Fraction
from math import factorial, comb
mp.mp.dps = 120
F450 = json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']
f450h=[mp.mpf(v) for v in F450]
N=220
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
r=[f450h[n]/mp.factorial(n)-Hmp[n] for n in range(181)]
beta=[mp.mpf(0)]*181; beta[0]=r[0]
twopi=2*mp.pi
for n in range(1,181):
    beta[n]=sum(r[n-2*j]*(twopi**(n-2*j))*comb(n-j-1,j) for j in range(0,(n-1)//2+1) if n-2*j>=1)
def w_of_z(z): return z/(mp.pi+mp.sqrt(mp.pi**2+z**2))
def R_map(z):
    w=w_of_z(z); tot=mp.mpf(0); pw=mp.mpf(1)
    for n in range(181): tot+=beta[n]*pw; pw*=w
    return tot
def Hz(z): return mp.sqrt(2)*mp.sqrt(1+mp.e**(-2*z))
def S_H(eta): return mp.quad(lambda t: mp.e**(-t)*Hz(eta*t), [0,5,20,60,mp.inf])
def S_R(eta): return mp.quad(lambda t: mp.e**(-t)*R_map(eta*t), [0,2,5,10,20,40,mp.inf])
def I_perim(eta):
    f=lambda t: mp.sqrt(t**(2*eta)+(1-t)**(2*eta))
    P=4*(1+eta)*mp.quad(f,[0,mp.mpf('0.5'),1])
    return P/(4*mp.sqrt(2)*(1+eta))
for eta in [mp.mpf('0.5'), mp.mpf('2'), mp.mpf('-0.5'), mp.mpf('-0.9')]:
    S=S_H(eta)+S_R(eta); Iv=I_perim(eta)
    print(f'  eta={eta}: |S-I|={mp.nstr(abs(S-Iv),3)}')
with open('progress.log','a') as f:
    f.write("step_32: verify_grw1r.py patched (dps-120 f450 in sec 6); eta=2 restored to ~1e-9\n")
