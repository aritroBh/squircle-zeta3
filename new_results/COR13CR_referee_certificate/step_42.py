import time, sys
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC2R_referee_certificate')
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
import mpmath as mp
mp.mp.dps = 65
import reflib as R
import cmzv6
T0 = time.monotonic()
# time 12 weight-4 words (cold cache)
ws = list(cmzv6.words(4))[:12]
for w in ws:
    R.cmzv_value(w)
el = (time.monotonic()-T0)/12
print(f"{el:.3f} s per weight-4 word (cold) -> full 1764 words ~ {el*1764/60:.1f} min")
