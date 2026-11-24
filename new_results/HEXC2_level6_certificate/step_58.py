import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
from hexc2_common import target4, atom_pi4, atom_ipi
from cmzv6 import mul
from fractions import Fraction as Fr
t0=time.monotonic()
T4 = target4()
print('T4 build: %.1fs, support %d'%(time.monotonic()-t0, len(T4)))
q0 = T4.get((1,1,1,1), Fr(0))
print('coefficient of word (1,1,1,1) in T4:', q0, '=', float(q0))
# closed form contributes -949/2160 * 1944 = ?
print('closed-form contribution alone:', Fr(-949,2160)*1944)
# pi^4 relation: (i pi)^4 = 81 * (1)^shuffle4 = 81*24*(1,1,1,1)
p = atom_ipi()
sh4 = mul(mul(p,p),mul(p,p))
print('(i pi)^4 as words:', sh4)
print('=> val((1,1,1,1)) = pi^4/1944; 1944*(coef) =', 1944*q0, 'should equal pi^4-part of T4 (numerically the value of T4 is 0)')
