code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, fsum
import json
setdps(80)
s = 4
def Jmod(u):
    # exact J via modular form (fast, exact to 1e-70 for u>0.4) 
    from mpmath import pi
    tot = pi**2*u/12 - log(2)/2 + 1/(24*u)
    n = 0
    while True:
        n += 1
        t = log(1+exp(-2*pi**2*u*n))
        tot -= t
        if t < mpf(10)**(-75): break
    return tot
def Jany(u):
    if u > 0.35: return Jmod(u)
    tot = mpf(0); d = 0
    while True:
        d += 1
        t = log(1+exp(-d/u))
        tot += t
        if t < mpf(10)**(-75): break
    return tot
p = 128
# Riemann sum (1/p) sum_k (k/p)^{-s} J(k/p)
rs = mpf(0)
K = 40*p
for k in range(1, K+1):
    u = mpf(k)/p
    rs += Jany(u)/u**s
rs /= p
print("Riemann sum (1/p) sum F(k/p):", nstr(rs, 25))
# Uq_direct value
data = json.load(open("Uq_data_s4.json"))
u128 = mpf(data["128"][0])
print("p^3 * Uq_direct(128)    :", nstr(u128*p**(s-1), 25))
print("continuum integral      : 2.2767747004318550048")
# also the actual inner(k) vs J(k/p) for p=128
def inner_direct(k, p):
    kf = mpf(k)
    return fsum([log(1+power(1-mpf(d)/kf, p)) for d in range(1, k)])
for k in [128, 256, 640]:
    print(f"k={k}: inner={nstr(inner_direct(k,p), 15)}  J(k/p)={nstr(Jany(mpf(k)/p), 15)}")
'''
exec(code)
