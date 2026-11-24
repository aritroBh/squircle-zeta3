
import mpmath as mp
mp.mp.dps = 50
sqrt2 = mp.sqrt(2); pi = mp.pi

def Hp3(t):
    # H'(t)/t^3 in manifestly regular form; H(t) = G(t) - sqrt2 t (rho-1)^2/rho^2 = lower-box margin /(N/2)
    q = mp.sqrt((1+sqrt2*t)/(1+sqrt2*t+t*t))
    w = 1/((1+sqrt2*t+t*t)*(1+q))
    Gp = (sqrt2/3)*(sqrt2+2*t)/(1+sqrt2*t+t*t)**2
    hp = sqrt2*w*w*t + sqrt2*w*t*(2+sqrt2*t)/(q*(1+sqrt2*t+t*t)**2)
    return Gp - hp

def Hfun(t):
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    psi = sqrt2*t - mp.log(1+sqrt2*t+t*t)
    G = psi - (sqrt2/3)*t**3/(1+sqrt2*t+t*t)
    return G - sqrt2*t*(rho-1)**2/rho**2

# verify H' = t^3 Hp3 numerically (finite difference) and positivity
mn = mp.mpf('9'); argmn = None
for k in range(20001):
    t = mp.mpf('0.127')*k/20000
    v = Hp3(t)
    if v < mn: mn, argmn = v, t
print("min Hp3 on [0,0.127] (20001 pts) =", mp.nstr(mn, 12), "at t =", argmn)
# check derivative relation
for tv in ['0.001','0.05','0.126']:
    t = mp.mpf(tv)
    fd = (Hfun(t+mp.mpf('1e-8')) - Hfun(t-mp.mpf('1e-8')))/mp.mpf('2e-8')
    print(f"t={tv}: H' fd = {mp.nstr(fd,10)}  t^3 Hp3 = {mp.nstr(t**3*Hp3(t),10)}")
# how far does Hp3 stay positive?
import numpy as np
last = None
for k in range(1, 10001):
    t = mp.mpf(k)*mp.mpf('0.001')
    if Hp3(t) < 0:
        last = t; break
print("Hp3 first negative near t =", last)
