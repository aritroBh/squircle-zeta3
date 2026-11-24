import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import mpmath as mp
import numeric_cmzv as NC
from hexc2_common import li_atom, cross_table
from fractions import Fraction as Fr
NC.set_dps(70)
def cmzv_combo(d):
    tot = mp.mpc(0)
    for w,c in d.items():
        tot += mp.mpf(c.numerator)/mp.mpf(c.denominator) * NC.cmzv_value(w)
    return tot
tab = cross_table()
need = {2: [Fr(1,3)], 3: [Fr(1,3), Fr(-1,2)], 4: [Fr(1,3), Fr(-1,3), Fr(2,3), Fr(1,2), Fr(-1,2), Fr(-2)]}
worst = 0; fails=[]
for n, zs in need.items():
    for z in zs:
        nreps = len(tab[z])
        conv = li_atom(n, z)
        got = cmzv_combo(conv)
        expect = mp.polylog(n, z)
        err = abs(got-expect)
        rel = err/max(1, abs(expect))
        worst = max(worst, float(rel))
        flag = 'OK' if rel < 1e-62 else 'FAIL'
        if rel >= 1e-62: fails.append((n,z,float(rel)))
        print(f'Li_{n}({z}) [{nreps} reps] err=%.2e %s'%(rel,flag))
print('worst:', worst, 'fails:', fails)
