import mpmath
from mpmath import mp, mpf, log, exp, quad, zeta, nstr
mp.dps = 100
A2 = quad(lambda t: log(1+t)**2/(2*t), [0, mpf('0.5'), 1])
B1 = quad(lambda s: log(1+exp(-2*s))**2, [0, 5, 20, 40])
print("A2 =", nstr(A2,80))
print("B1 =", nstr(B1,80))
print("B1-A2 =", nstr(B1-A2,6))
# manual midpoint sum with tiny step on [0,30], tail < e^{-60}
mp.dps = 50
h = mpf('0.00025')
tot = mpf(0)
N = int(30/h)
s = h/2
for i in range(N):
    tot += log(1+exp(-2*s))**2
    s += h
tot *= h
print("manual:", nstr(tot,50))
print("zeta3/8:", nstr(zeta(3)/8,50))
print("manual - zeta3/8:", nstr(tot-zeta(3)/8,6))
