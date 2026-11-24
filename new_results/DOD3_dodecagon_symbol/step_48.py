import sys, os
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 620
# L(3, chi_12) via Hurwitz zeta (exact evaluation)
L3c = (mp.zeta(3, mp.mpf(1)/12) - mp.zeta(3, mp.mpf(5)/12) - mp.zeta(3, mp.mpf(7)/12) + mp.zeta(3, mp.mpf(11)/12)) / 12**3
print("L(3,chi_12) =", mp.nstr(L3c, 40))
# verify vs partial series
mp.mp.dps = 100
chi = {1:1, 5:-1, 7:-1, 11:1}
s = mp.mpf(0)
for n in range(1, 400001):
    c = chi.get(n%12, 0)
    if c: s += c/mp.power(n,3)
print("series check |diff|:", mp.nstr(abs(s-L3c),3))
mp.mp.dps = 620
val = mp.mpf('-12.00785834784103590809332')  # need full 620-digit value; recompute below properly
