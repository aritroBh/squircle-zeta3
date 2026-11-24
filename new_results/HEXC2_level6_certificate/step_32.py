import sys; sys.dont_write_bytecode = True; sys.path.insert(0,'.')
from hexc2_common import li_atom, cross_table, FRNAM, MAPS
from fractions import Fraction as Fr
import math
val = {1: 1j*math.pi/3, 2: -math.log(3)/2 + 1j*math.pi/6, 3: -math.log(2),
       4: -math.log(3)/2 - 1j*math.pi/6, 5: -1j*math.pi/3}
def ev(d):
    tot = 0
    for w,c in d.items():
        assert len(w)==1
        tot += float(c)*val[w[0]]
    return tot
tab = cross_table()
z = Fr(-1,3)
for i, rep in enumerate(tab[z]):
    e,p0,p1,pinf = rep
    EXP = {Fr(0):0, Fr(1,2):1, Fr(1):2, Fr(-1):5, None:4}
    print('rep', i, tuple(FRNAM[q] for q in rep), ' u-path exps:', EXP[p0], '->', EXP[e], ' sing exps:', EXP[p1], EXP[pinf])
    print('   value:', ev(li_atom(1, z, rep=i)), ' expect', math.log(4/3))
