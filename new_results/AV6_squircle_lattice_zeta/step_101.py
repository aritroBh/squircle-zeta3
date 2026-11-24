code = r'''
from c2 import *
from squircle import setdps
from mpmath import mp, nstr, log, power, fsum, mpf
setdps(60)
s = 4
# brute G_l(4) to K=1200
P = mpf(0)
acc = mpf(0)
for k in range(1, 1201):
    acc += k*log(k)
    if k >= 2: P += power(k, -s-1)*log(k)*acc
print("brute G_l(4) K=1200 (no tail):", nstr(P, 20))
# my G_l(4) with tail:
print("G_l(4) full:", nstr(G_l(s), 20))
# compute the tail alone
tail = G_l(s) - P
print("tail computed:", nstr(tail, 12))
# estimate tail independently: sum_{m<=1200} m ln m * HTln(1200, 5)
from c2 import HTln
est = HTln(1200, s+1)*fsum([m*log(m) for m in range(1, 1201)])
print("est first part:", nstr(est, 12))
'''
exec(code)
