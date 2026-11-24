"""w6 identification step 2a: expensive constants, computed once at dps 270 and cached as strings in consts_270.json.
  T  = sum_{n>=1} r_n / n^5,  r_n = sum_{m>=n} (-1)^{m-n}/m = int_0^1 x^(n-1)/(1+x) dx
     = zeta(6) - int_0^1 Li5(x)/(1+x) dx
  S6 = zeta(-5,-1) = sum_{n>m>=1} (-1)^(n+m)/(n^5 m) = log2 * eta(5) - T,   eta(5) = (15/16) zeta(5)
Low-precision self-check of S6 against the defining double sum is done in s09b (printed there).
"""
import sys, os, json, time
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
fn = os.path.join(W, "consts_270.json")
if os.path.exists(fn):
    print("exists"); sys.exit(0)
mp.mp.dps = 280          # after imports
print("mp.mp.dps =", mp.mp.dps)
t0 = time.time()
I = mp.quad(lambda x: mp.polylog(5, x) / (1 + x), [0, mp.mpf(1) / 2, 1], error=True, maxdegree=12)
print("quad I =", mp.nstr(I[0], 40), "err", mp.nstr(I[1], 3), f"[{time.time()-t0:.0f}s]")
T = mp.zeta(6) - I[0]
S6 = mp.log(2) * mp.mpf(15) / 16 * mp.zeta(5) - T
json.dump({"dps": 280, "T": mp.nstr(T, 275), "S6": mp.nstr(S6, 275), "quad_err": mp.nstr(I[1], 5)}, open(fn, "w"), indent=1)
print("S6 =", mp.nstr(S6, 50))
