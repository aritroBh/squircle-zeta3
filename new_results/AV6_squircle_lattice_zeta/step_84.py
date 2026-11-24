code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, fsum, pi, quad, inf
setdps(80)
s = 4
def M2(u): return fsum([d*d/(exp(d/u)+1) for d in range(1, 80)])
def A1(u): return -M2(u)/(2*u*u)
def F1(u): return A1(u)/u**s
# integral of F1
I_F1 = -gamma(s+1)*(1-power(2,-s))*zeta(s+1)*zeta(s-1)/2
print("int F1 =", nstr(I_F1, 20))
# Riemann sums at several p with big cutoff
for p in [128, 256, 512]:
    rs = mpf(0)
    K = 60*p
    for k in range(1, K+1):
        rs += F1(mpf(k)/p)
    rs /= p
    # tail integral estimate beyond u=60: A1(u) ~ -0.9 u -> F1 ~ -0.9/u^3 : int_60^inf ~ -0.9/(2*3600)
    print(f"p={p}: Riem={nstr(rs, 20)}  err={(nstr(rs - I_F1, 6))}  err*p={nstr((rs-I_F1)*p, 8)}")
'''
exec(code)
