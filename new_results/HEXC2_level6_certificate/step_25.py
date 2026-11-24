import sys; sys.dont_write_bytecode = True; sys.path.insert(0,'.')
from hexc2_common import li_atom
from fractions import Fraction as Fr
# numeric spot check of letter values
# (j) = -log(1-zeta^j): (1)=i pi/3, (2)=-log3/2 + i pi/6, (3)=-log2, (4)=-log3/2 - i pi/6, (5)=-i pi/3
import cmath, math
val = {1: 1j*math.pi/3, 2: -math.log(3)/2 + 1j*math.pi/6, 3: -math.log(2),
       4: -math.log(3)/2 - 1j*math.pi/6, 5: -1j*math.pi/3}
def ev(d):
    tot = 0
    for w,c in d.items():
        assert len(w)==1
        tot += float(c)*val[w[0]]
    return tot
print('Li_1(1/3) should be', math.log(1.5), '; got', ev(li_atom(1, Fr(1,3))))
print('Li_1(1/2) should be', math.log(2), '; got', ev(li_atom(1, Fr(1,2))))
print('Li_1(-1/2) should be', math.log(2/3), '; got', ev(li_atom(1, Fr(-1,2))))
print('Li_1(-2) should be', -math.log(3), '; got', ev(li_atom(1, Fr(-2))))
print('Li_1(2/3) should be', math.log(3), '; got', ev(li_atom(1, Fr(2,3))))
