"""w6 identification step 2b: low-precision check that S6 (cached) equals the defining alternating double sum
sum_{n>m>=1} (-1)^(n+m)/(n^5 m), by direct partial summation to N terms plus the explicit tail bound.
Also checks S6 at a second precision (recompute with dps 120 quad) for the cache's reliability."""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
mp.mp.dps = 40
c = json.load(open(os.path.join(W, "consts_270.json")))
S6 = mp.mpf(c["S6"])
N = 200000
tot = mp.mpf(0); a = mp.mpf(0)      # a = sum_{m<n} (-1)^m/m
# sum_{n=2}^N (-1)^n a_{n-1}/n^5
a = -mp.mpf(1)
sgn = 1
tot = mp.mpf(0)
for n in range(2, N + 1):
    tot += (-1) ** n * a / mp.mpf(n) ** 5
    a += (-1) ** n / mp.mpf(n)
print("partial sum N =", N, ":", mp.nstr(tot, 25))
print("cached S6          :", mp.nstr(S6, 25))
print("difference:", mp.nstr(tot - S6, 5), " (tail bound ~ 1/(4 N^4) =", mp.nstr(mp.mpf(1) / (4 * mp.mpf(N) ** 4), 3), ")  dps", mp.mp.dps)
# second-precision recompute
mp.mp.dps = 130
I = mp.quad(lambda x: mp.polylog(5, x) / (1 + x), [0, mp.mpf(1) / 2, 1])
S6b = mp.log(2) * mp.mpf(15) / 16 * mp.zeta(5) - (mp.zeta(6) - I)
print("S6 (dps 130 recompute) - cached:", mp.nstr(S6b - mp.mpf(c["S6"]), 5), " dps", mp.mp.dps)
