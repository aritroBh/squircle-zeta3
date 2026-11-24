
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 65
import reflib as R
t0 = time.monotonic()
vals = [R.cmzv_value(w) for w in R.words(4)[:20]]
t1 = time.monotonic()
print(f"20 weight-4 words: {t1-t0:.2f}s  -> {(t1-t0)/20*1000:.0f} ms/word")
t0 = time.monotonic()
vals = [R.cmzv_value(w) for w in R.words(3)]
t1 = time.monotonic()
print(f"all 252 weight-3 words: {t1-t0:.2f}s")
# spot re-check after optimization
print("zeta2 check:", mp.nstr(R.cmzv_value((-1,0)) - mp.zeta(2), 3))
print("w4 sample:", mp.nstr(vals[0], 10))
