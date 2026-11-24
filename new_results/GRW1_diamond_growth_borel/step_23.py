
import mpmath as mp
mp.mp.dps = 90
import json
from fractions import Fraction
from math import factorial as ffact

N = 260
KK = N//2+2
a = [Fraction(1, ffact(2*k)) for k in range(KK)]
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hf = [Fraction(0)]*(N+1)
for n in range(N+1):
    Hf[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*ffact(n-2*k)) for k in range(n//2+1))
H = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hf]
Hregf = {}
def Hreg(j, x):
    tot = mp.mpf(0); mx = mp.mpf(0)
    for m in range(N-j+1):
        t = H[m]*x**(m+j)/mp.factorial(m+j)
        tot += t; mx = max(mx, abs(t))
        if m > 2*j+20 and abs(t) < mp.mpf('1e-60')*mx and abs(t) < abs(tot)*mp.mpf('1e-60'):
            break
    return tot

def kj(j, u):
    L = -mp.log(1+mp.exp(-u))
    return (1+mp.exp(-u))**(-2) * L**j / mp.factorial(j)

def de_quad(phi, h, Tlo, Thi):
    tot = mp.mpf(0)
    klo = int(Tlo/h)+1; khi = int(Thi/h)+1
    for k in range(-klo, khi+1):
        t = k*h
        sh = mp.sinh(t); u = mp.exp(mp.pi/2*sh)
        up = mp.pi/2*mp.cosh(t)*u
        tot += phi(u)*up
    return tot*h

Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e20')))  # u from 1e-20 to 36
Thi = mp.asinh(2/mp.pi*mp.log(36))
def BBj(j, z):
    if j == 0:
        phi = lambda u: mp.exp(-u)*kj(0,u)*Hreg(0, z*u)
    else:
        phi = lambda u: mp.exp(-u)*kj(j,u)*u**(-j)*Hreg(j, z*u)
    return de_quad(phi, mp.mpf(1)/40, Tlo, Thi)

def B_exact(z, jmax=12):
    return sum(BBj(j, z) for j in range(0, jmax+1))

P = json.load(open("pade9090.json"))
num = [mp.mpf(v) for v in P['num']]; den = [mp.mpf(v) for v in P['den']]
def peval(t):
    nv = num[0]; dv = den[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv

d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
def B_series(z):
    return sum(f[n]/mp.factorial(n)*z**n for n in range(181))

print("=== B(zeta) real axis: exact Bj-sum vs Pade vs series ===")
for z in [mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('3.0'), mp.mpf('-1.0'), mp.mpf('-3.0')]:
    be = B_exact(z)
    bp = peval(z)
    line = f"z={mp.nstr(z,6)}: Bexact={mp.nstr(be,25)}"
    print(line)
    line2 = f"          Pade ={mp.nstr(bp,25)}  |diff|={mp.nstr(abs(be-bp),3)}"
    if abs(z) < mp.pi/2:
        bs = B_series(z)
        line2 += f"  series={mp.nstr(bs,20)} |d|={mp.nstr(abs(be-bs),3)}"
    print(line2, flush=True)
