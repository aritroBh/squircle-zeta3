
src = open("reflib.py").read()
src = src.replace("i1 = [inv[n][s[0] - 1] for n in range(K + 1)]",
                  "i1 = [mp.mpf(1)] + [inv[n][s[0] - 1] for n in range(1, K + 1)]")
src = src.replace("i2 = [inv[n][s[1] - 1] for n in range(K + 1)]",
                  "i2 = [mp.mpf(1)] + [inv[n][s[1] - 1] for n in range(1, K + 1)]")
src = src.replace("i3 = [inv[n][s[2] - 1] for n in range(K + 1)]",
                  "i3 = [mp.mpf(1)] + [inv[n][s[2] - 1] for n in range(1, K + 1)]")
src = src.replace("i4 = [inv[n][s[3] - 1] for n in range(K + 1)]",
                  "i4 = [mp.mpf(1)] + [inv[n][s[3] - 1] for n in range(1, K + 1)]")
open("reflib.py", "w").write(src)
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 65
import reflib as R
t0 = time.monotonic()
vals = [R.cmzv_value(w) for w in R.words(4)[:20]]
t1 = time.monotonic()
print(f"20 weight-4 words: {t1-t0:.2f}s")
t0 = time.monotonic()
vals3 = [R.cmzv_value(w) for w in R.words(3)]
t1 = time.monotonic()
print(f"252 weight-3 words: {t1-t0:.2f}s")
print("zeta2 check:", mp.nstr(R.cmzv_value((-1,0)) - mp.zeta(2), 3))
