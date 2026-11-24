
import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))

# EXACT regular form: Btrue = rho*Q/(t^2 (2+sqrt2 t)),
#   Q = sqrt2 t (1-rho) + t^2/2 + (1-rho)^2 - rho^2 A1', A1' = 1 - atan(zz)/zz, zz = sqrt2 t/(2+sqrt2 t)
def Btrue_form(t):
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    zz = sqrt2*t/(2+sqrt2*t)
    Ap = 1 - mp.atan(zz)/zz if t != 0 else mp.mpf('0')
    Q = sqrt2*t*(1-rho) + t*t/2 + (1-rho)**2 - rho**2*Ap
    return rho*Q/(t*t*(2+sqrt2*t))

worst = mp.mpf('0')
import random
random.seed(7)
for _ in range(600):
    t = mp.mpf(10)**random.uniform(-12, 0)   # t in [1e-12, 1], covers final cell and [100,400] range
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    phi = mp.im(Phi(x0,N)) - th
    Bt = (phi + pi/4)/((pi/sqrt2)*y)
    Bf = Btrue_form(t)
    worst = max(worst, abs(Bt-Bf)/abs(Bt))
print("EXACT FORM Btrue = rho Q/(t^2(2+rt2 t)): worst rel resid over t in [1e-12,1] =", mp.nstr(worst, 8))

# sup |Btrue| on the final cell t in [0, T], T = (1+y_J^2)*y_J, y_J = sqrt(pi/2 e-10)
yJ = mp.sqrt(pi/2*mp.mpf('1e-10')); T = (1+yJ**2)*yJ
print("y_J =", mp.nstr(yJ, 12), " T =", mp.nstr(T,12))
import numpy as np
mx = mp.mpf('0')
for k in range(20001):
    t = T*k/20000
    v = abs(Btrue_form(t)) if t>0 else mp.mpf(1)/6
    mx = max(mx, v)
print("sup |Btrue| on [0,T] (dense grid, 20001 pts) =", mp.nstr(mx, 15), " 1/6 =", mp.nstr(mp.mpf(1)/6,15))

# also sup on the [100,400] domain t range: N in [100,400] => y in [sqrt(pi/800), sqrt(pi/200)]
y1 = mp.sqrt(pi/800); y2 = mp.sqrt(pi/200)
t1 = y1; t2 = y2*(1+y2**2)
mx2 = mp.mpf('0')
for k in range(20001):
    t = t1 + (t2-t1)*k/20000
    mx2 = max(mx2, abs(Btrue_form(t)))
print("sup |Btrue| on [100,400]-domain t-range =", mp.nstr(mx2, 15))
