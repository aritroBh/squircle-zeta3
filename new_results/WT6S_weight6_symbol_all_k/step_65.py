
import flint
from fractions import Fraction as F
M = flint.fmpq_mat([['0','0','1','1/4'],['0','3/16','0','-1/8'],['0','0','0','1/8'],['3/4','0','-3','3/4']])
d = M.det()
print("det:", d, type(d))
print("d == F(9,512):", d == F(9,512))
print("F(str(d)) == F(9,512):", F(str(d)) == F(9,512))
print("d == flint.fmpq('9/512'):", d == flint.fmpq('9/512'))
