from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, im, re, sqrt, zeta,
                    gammainc, hyper, quad, log, nsum, inf, nprint, fac)
mp.dps = 60

def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

# ---------- Fourier coefficients of f_s = |x|^s (cosine) and g_s = sgn(x)|x|^s (sine) on [-pi,pi] ----------
def a_cos(n, s):   # (2/pi) int_0^pi x^s cos(nx) dx
    return (2/pi)*quad(lambda x: x**s*cos(n*x), [0, pi])
def b_sin(n, s):   # (2/pi) int_0^pi x^s sin(nx) dx
    return (2/pi)*quad(lambda x: x**s*sin(n*x), [0, pi])

# Closed form 1: hypergeometric 1F2
def a_cos_1F2(n, s):
    return (2/pi)* pi**(s+1)/(s+1) * hyper([ (s+1)/2 ], [ mpf('0.5'), (s+3)/2 ], -(n*pi)**2/4)
def b_sin_1F2(n, s):
    return (2/pi)* n*pi**(s+2)/(s+2) * hyper([ (s+2)/2 ], [ mpf('1.5'), (s+4)/2 ], -(n*pi)**2/4)

# Closed form 2: lower incomplete gamma  int_0^pi x^s e^{inx} dx = (-in)^(-s-1) gamma(s+1, -i n pi)
def Iexp(n, s):
    return (-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi)
def a_cos_gam(n, s): return re((2/pi)*Iexp(n,s))
def b_sin_gam(n, s): return im((2/pi)*Iexp(n,s))

print("== coefficient closed-form cross-checks (60 dps) ==")
for s in [mpf('2.5'), mpf('3'), mpf('3.7')]:
    for n in [1,2,5,17]:
        v0 = a_cos(n,s); v1 = a_cos_1F2(n,s); v2 = a_cos_gam(n,s)
        w0 = b_sin(n,s); w1 = b_sin_1F2(n,s); w2 = b_sin_gam(n,s)
        print(f"s={s}, n={n}: a err1F2={abs(v0-v1):.2e} errGam={abs(v0-v2):.2e} | b err1F2={abs(w0-w1):.2e} errGam={abs(w0-w2):.2e}")
logline("step_04: coefficient closed forms 1F2 & incomplete-gamma verified vs quadrature, s=2.5,3,3.7, n<=17, 60dps")
