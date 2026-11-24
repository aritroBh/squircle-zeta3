import time, sys, json
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/COR13C_third_corner_certificate')
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
import mpmath as mp
from fractions import Fraction as Fr
import reflib13 as R13
import reflib as R
import cor13c_common as HC
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# compare frame choices for each Li argument
for z in ("1/3", "-1/3", "2/3", "1/2", "-1/2", "-2"):
    hc = HC.cross_table()[Fr(z)]
    mine = R13.cross_reps13(Fr(z))
    log(f"z={z}: HC rep0 frame {hc[0]} vs my rep0 frame {mine[0]}; same: {hc[0] == mine[0]}; #frames {len(hc)}/{len(mine)}")
