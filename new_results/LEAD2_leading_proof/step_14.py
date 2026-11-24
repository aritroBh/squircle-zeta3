import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, pi, zeta, gamma, sin, quad
import sys; sys.path.insert(0, '.')
from lead2_lib import E_coeffs, Efun, m_e, m_odd, L1, L2

E = E_coeffs(20)
lam = log(2)
# d=1 quantities
c2 = lam**2/4 + pi**2/24
K1 = c2  # delta_1 = 1
e2 = E[2]
print("K1 =", mpmath.nstr(K1,30), " target (1!/2)e2 =", mpmath.nstr(e2/2,30))

def A(p,q):
    return (-1)**p*quad(lambda r: r**p*log(1-mpmath.exp(-r))**q/(mpmath.exp(r)-1), [0,1,mpmath.inf])
T1 = -A(1,0) + mpf('0.5')*A(1,0)   # = -(1/2)A_{1,0}
print("T1 =", mpmath.nstr(T1,30), " pi^2/12 =", mpmath.nstr(pi**2/12,30))
piece2_1 = -mpf('0.5')*quad(lambda u: (L1(u)-L2(u))*u*m_e(u), [-1,0,1])
print("piece2(1) direct =", mpmath.nstr(piece2_1,30))
print("T1+piece2 =", mpmath.nstr(T1+piece2_1,30), "vs K1")
print("pi^2/3 - 2e2 =", mpmath.nstr(pi**2/3-2*e2,30), " /4 =", mpmath.nstr((pi**2/3-2*e2)/4,30))

# strip GFs
def hfun(w): return ((1+w)/sqrt(1+w**2)-1)/w
def qfun(w): return (1+w)/(1+w**2)**mpf('1.5')
def g1fun(w): return (1-w)/(1+w)*hfun(w)
for s in [mpf('0.3'), mpf('-0.4'), mpf('0.73')]:
    Fcal = quad(lambda w: w**(-s)*hfun(w), [0,1,mpmath.inf])
    rhsF = (Efun(s)-Efun(-s))/s
    Wcal = quad(lambda w: w**(-s)*qfun(w), [0,1,mpmath.inf])
    rhsW = Efun(s)+Efun(-s)
    Gcal = quad(lambda w: w**(-s)*g1fun(w), [0,1,mpmath.inf])
    rhsG = 2*pi/sin(pi*s) - (Efun(s)+Efun(-s))/s
    print(f"s={s}: F: {mpmath.nstr(Fcal-rhsF,4)}   W: {mpmath.nstr(Wcal-rhsW,4)}   G: {mpmath.nstr(Gcal-rhsG,4)}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: d=1 decomposition + strip GF identities (F,E,W,G) numerics\n")
