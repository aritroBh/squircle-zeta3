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
    return fsum([d*d/(exp(d/u)+1) for d in range(1, 60)])
def S3(u):
    return fsum([d**3/(exp(d/u)+1) for d in range(1, 60)])
def T4(u):
    return fsum([d**4/(exp(d/u)+1)**2 for d in range(1, 60)])
p = 128
# Riemann sums
rsJ = rsA1 = rsA2 = mpf(0)
K = 40*p
for k in range(1, K+1):
    u = mpf(k)/p
    w = u**(-s)/p
    rsJ += w*J(u)
    rsA1 += w*(-M2(u)/(2*u*u))
    rsA2 += w*(-S3(u)/(3*u**3) + T4(u)/(8*u**4))
print("Riem J          :", nstr(rsJ, 20))
print("Riem A1/p       :", nstr(rsA1/p, 20))
print("Riem A2/p^2     :", nstr(rsA2/p**2, 20))
print("total           :", nstr(rsJ + rsA1/p + rsA2/p**2, 20))
print("data D(128)     : 2.182068885825243373674029")
# A1 integral value
print("A1 integral/128 :", nstr(-gamma(s+1)*eta(s+1)*zeta(s-1)/(2*p), 20))
'''
exec(code)
