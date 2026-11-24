import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import mpmath as mp
mp.mp.dps = 150
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
    return tot, mx
for x in [-80, -200]:
    v, mx = Hb0(mp.mpf(x))
    print(f"H^[0]({x}) = {mp.nstr(v, 20)}   peak term {mp.nstr(mx, 4)}   2e^(2sqrt|x|) = {mp.nstr(2*mp.e**(2*mp.sqrt(-x)), 6)}")

# direct keyholed-circle contour integral at x = -200: circle |s| = sqrt(X) plus loops.
# Use the PLAIN circle at a radius between consecutive |s_j|: choose R = sqrt(X) adjusted to avoid |s_j|
X = mp.mpf(200)
R = mp.sqrt(X)
# |s_j| = 2X/(pi|2k+1|): find which are near R
print(" |s_j| values near R =", mp.nstr(R, 8), ":", [mp.nstr(2*X/(mp.pi*abs(2*k+1)), 6) for k in range(-6, 6)])
# pick R mid-way between two consecutive |s_j|
R = (2*X/(mp.pi*9) + 2*X/(mp.pi*11))/2
print("using R =", mp.nstr(R, 8))
x = -X
iu = mp.mpc(0,1)
Hfun = lambda z: mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))
def integrand_phi(phi):
    s = R*mp.e**(iu*phi)
    return mp.e**s * Hfun(x/s)/s * (iu*R*mp.e**(iu*phi))
cont = mp.quad(integrand_phi, [0, mp.pi/2, mp.pi, 3*mp.pi/2, 2*mp.pi])/(2*mp.pi*iu)
print("plain circle (no keyholes) contour =", mp.nstr(cont, 18))
'''
with open(wd+"/step_37.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_37.py"], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-2500:]); print("STDERR:", r.stderr[-600:]); print("elapsed", time.time()-t0)
