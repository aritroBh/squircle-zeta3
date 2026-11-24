code = r'''
from c2 import *
from squircle import setdps
from mpmath import mp, nstr, log, power, fsum, mpf
setdps(60)
print("G_2(5) =", nstr(G_2(5), 25))
print("G_x(5) =", nstr(G_x(5), 25))
print("sum    =", nstr(G_2(5)+G_x(5), 25), " (expect 0.0321016)")
# brute G_2(5)
P = mpf(0); acc = mpf(0)
for k in range(1, 1501):
    if k >= 2: P += power(k, -6)*acc
    acc += k*log(k)**2
print("G_2(5) brute partial:", nstr(P, 20))
# brute G_x(5)
Px = mpf(0)
for k in range(2, 801):
    sx = fsum([m*log(m)*(k-m)*log(k-m) for m in range(1, k)])
    Px += power(k, -6)*sx
print("G_x(5) brute partial:", nstr(Px, 20))
'''
exec(code)
