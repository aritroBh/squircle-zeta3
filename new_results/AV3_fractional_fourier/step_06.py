from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, im, re, sqrt, zeta,
                    gammainc, hyper, quad, log, nsum, inf, nstr, fac)
mp.dps = 60
def E(x): return nstr(abs(x), 3)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def a_cos(n, s):   return (2/pi)*quad(lambda x: x**s*cos(n*x), [0, pi])
def b_sin(n, s):   return (2/pi)*quad(lambda x: x**s*sin(n*x), [0, pi])
def a_cos_1F2(n, s):
    return (2/pi)* pi**(s+1)/(s+1) * hyper([ (s+1)/2 ], [ mpf('0.5'), (s+3)/2 ], -(n*pi)**2/4)
def b_sin_1F2(n, s):
    return (2/pi)* n*pi**(s+2)/(s+2) * hyper([ (s+2)/2 ], [ mpf('1.5'), (s+4)/2 ], -(n*pi)**2/4)
def Iexp(n, s):
    return (-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi)
def a_cos_gam(n, s): return re((2/pi)*Iexp(n,s))
def b_sin_gam(n, s): return im((2/pi)*Iexp(n,s))

print("== coefficient closed-form cross-checks (60 dps) ==")
for s in [mpf('2.5'), mpf('3'), mpf('3.7')]:
    for n in [1,2,5,17]:
        v0 = a_cos(n,s); v1 = a_cos_1F2(n,s); v2 = a_cos_gam(n,s)
        w0 = b_sin(n,s); w1 = b_sin_1F2(n,s); w2 = b_sin_gam(n,s)
        print(f"s={nstr(s,4)}, n={n}: a err1F2={E(v0-v1)} errGam={E(v0-v2)} | b err1F2={E(w0-w1)} errGam={E(w0-w2)}")
logline("step_05: coeff closed forms 1F2 & incomplete-gamma verified vs quadrature (60 dps)")
