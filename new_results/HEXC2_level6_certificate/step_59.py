import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
t0 = time.monotonic()
from hexc2_common import target2, target3, target4, hpl_to_cmzv, li_atom
T2, T3, T4 = target2(), target3(), target4()
print('fresh target builds: %.1fs'%(time.monotonic()-t0))
import mpmath as mp
import numeric_cmzv as NC
NC.set_dps(62)
t1 = time.monotonic()
def combo(d):
    tot = mp.mpc(0)
    for w,c in d.items():
        tot += mp.mpf(c.numerator)/c.denominator * NC.cmzv_value(w)
    return tot
v4 = combo(T4)
print('T4 numeric at 62dps: %.1fs, value=%s'%(time.monotonic()-t1, mp.nstr(v4,6)))
