import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import mpmath as mp
mp.mp.dps = 100
from fractions import Fraction
from math import factorial as ffact

N = 420
KK = N//2+2
a = [Fraction(1, ffact(2*k)) for k in range(KK)]
s_ = [Fraction(0)]*KK; s_[0] = Fraction(1)
for k in range(1, KK):
    s_[k] = (a[k] - sum(s_[i]*s_[k-i] for i in range(1, k)))/2
H = [Fraction(0)]*(N+1)
for n in range(N+1):
    H[n] = 2*sum(s_[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*ffact(n-2*k)) for k in range(n//2+1))
def Hb0(x):
    tot = mp.mpf(0); mx = mp.mpf(0)
    for m in range(N+1):
        t = mp.mpf(H[m].numerator)/mp.mpf(H[m].denominator)*x**m/mp.factorial(m)
        tot += t; mx = max(mx, abs(t))
        if m > 60 and abs(t) < mp.mpf('1e-80')*mx: break
    return tot

# J_0-form for x <= 0: H^[0](x) = sqrt2 sum_{k>=0} c_k J_0(2 sqrt((1-2k) x))
ck = [mp.mpf(1)]
for k in range(1, 40001):
    ck.append(ck[-1]*(mp.mpf('0.5')-(k-1))/k)
def Jform(x, kmax=40000):
    tot = mp.mpf(0)
    for k in range(kmax+1):
        tot += ck[k]*mp.besselj(0, 2*mp.sqrt((1-2*k)*x))
    return mp.sqrt(2)*tot

print("x, H^[0](x) [series], J-form, diff, sqrt2*(I0(2sqrt|x|)+1) bound")
for x in [-1, -5, -20, -80, -200]:
    v1 = Hb0(mp.mpf(x))
    v2 = Jform(mp.mpf(x))
    bd = mp.sqrt(2)*(mp.besseli(0, 2*mp.sqrt(-x))+1)
    print(f"x={x}: series={mp.nstr(v1,16)}")
    print(f"       Jform ={mp.nstr(v2,16)}  |diff|={mp.nstr(abs(v1-v2),3)}  bound={mp.nstr(bd,6)}  ok={abs(v1)<=bd}")
'''
with open(wd+"/step_38.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_38.py"], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-2500:]); print("STDERR:", r.stderr[-600:]); print("elapsed", time.time()-t0)
