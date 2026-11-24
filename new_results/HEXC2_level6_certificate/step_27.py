import sys; sys.dont_write_bytecode = True; sys.path.insert(0,'.')
import importlib
import hexc2_common
importlib.reload(hexc2_common)
from hexc2_common import li_atom
from fractions import Fraction as Fr
import math
val = {1: 1j*math.pi/3, 2: -math.log(3)/2 + 1j*math.pi/6, 3: -math.log(2),
       4: -math.log(3)/2 - 1j*math.pi/6, 5: -1j*math.pi/3}
def ev(d):
    tot = 0
    for w,c in d.items():
        assert len(w)==1, w
        tot += float(c)*val[w[0]]
    return tot
for z, expect in [(Fr(1,3), math.log(1.5)), (Fr(1,2), math.log(2)), (Fr(-1,2), math.log(2/3)),
                  (Fr(-2), -math.log(3)), (Fr(2,3), math.log(3)), (Fr(-1,3), math.log(4/3))]:
    got = ev(li_atom(1, z))
    print(f'Li_1({z}) expect {expect:+.12f} got {got:+.12f}  {"OK" if abs(got-expect)<1e-12 else "FAIL"}')
