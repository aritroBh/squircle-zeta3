
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
s = open(os.path.join(WD, 'verify_wt6s.py')).read()
s = s.replace("    okD &= M.det() == F(9, 512)", "    okD &= F(str(M.det())) == F(9, 512)")
s = s.replace("""    combo = [sum((Q(-12, 209) * avec['Li_{3,3}(b,1)'][t],), ZEROQ) for t in range(1)]  # placeholder
    # M5b N6 control""", "    # M5b N6 control")
open(os.path.join(WD, 'verify_wt6s.py'), 'w').write(s)
# quick standalone test of just D logic
import subprocess
PY = '<project>/.venv/bin/python'
test = r'''
import json, pickle
from fractions import Fraction as F
import flint, qsqrt2
from w4_symbol import Q, R, B, ONE as ONEQ, ban_float, self_check
import wt6s_core as C6
from wt6s_core import depth2_symbol, pairing_vector6, add, scale, ZEROQ
ban_float(); self_check()
ATOMS = [("Li33b1", 3, 3, B, ONEQ), ("Li33s1", 3, 3, R/2, ONEQ),
         ("Li33b2", 3, 3, B, Q(2)), ("Li42bs", 4, 2, B, R/2)]
avec = {}
for name, m, n, x, y in ATOMS:
    avec[name] = pairing_vector6(depth2_symbol(m, n, x, y))
M = flint.fmpq_mat([[str(avec[nm][t].p) for nm, *_ in ATOMS] for t in range(4)])
print("det ok:", F(str(M.det())) == F(9, 512))
N6 = {"Li33b1": Q(-12, 209), "Li33s1": Q(-160, 80), "Li33b2": Q(0, F(329, 4)), "Li42bs": Q(0, 120)}
c5 = [sum((N6[nm]*avec[nm][t] for nm, *_ in ATOMS), ZEROQ) for t in range(4)]
print("N6 control:", c5 == [Q(0, F(449,4)), Q(-30), Q(0, 15), Q(-9)])
'''
open(os.path.join(WD, '_testD.py'), 'w').write(test)
r = subprocess.run([PY, '-B', '_testD.py'], capture_output=True, text=True, cwd=WD)
print(r.stdout, r.stderr[-500:])
os.remove(os.path.join(WD, '_testD.py'))
