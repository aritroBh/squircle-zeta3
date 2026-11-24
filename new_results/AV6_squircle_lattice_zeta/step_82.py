code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, fsum, pi
setdps(80)
s = 4
def J(u):
    tot = pi**2*u/12 - log(2)/2 + 1/(24*u)
    n = 0
    while True:
        n += 1
        t = log(1+exp(-2*pi**2*u*n))
        tot -= t
        if t < mpf(10)**(-75): break
    return tot
def M2(u):
    return fsum([d*d/(exp(d/u)+1) for d in range(1, 80)])
def inner(k, p):
    kf = mpf(k)
    return fsum([log(1+power(1-mpf(d)/kf, p)) for d in range(1, k)])
p = 256
print("u        inner-J         (inner-J)*p      -M2/(2u^2)      ratio")
for k in [64, 128, 256, 384, 512, 768, 1024]:
    u = mpf(k)/p
    d0 = inner(k, p) - J(u)
    a1 = -M2(u)/(2*u*u)
    print(f"u={nstr(u,5)}: {nstr(d0, 12)}  {nstr(d0*p, 12)}   {nstr(a1, 12)}   {nstr(d0*p/a1 if a1 else 0, 8)}")
'''
exec(code)
