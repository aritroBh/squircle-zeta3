import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import mpmath as mp
import fast_cmzv as FC
from hexc2_common import target2, target3, target4
FC.setup(62); FC.finish_setup()
mp.mp.dps = 65
def combo_fast(d):
    tot = mp.mpc(0)
    for w,c in d.items():
        v = FC.cmzv_value(w)
        tot += mp.mpf(c.numerator)/c.denominator * mp.mpc(v.real.mid(), v.imag.mid())
    return tot
for name, T in [('T2', target2), ('T3', target3), ('T4', target4)]:
    t0 = time.monotonic()
    v = combo_fast(T())
    print(f'{name}: {time.monotonic()-t0:.1f}s  |value| = {mp.nstr(abs(v),3)}')
