
import flint
from fractions import Fraction as F
M = flint.fmpq_mat([[0,0,1,F(1,4)],[0,F(3,16),0,F(-1,8)],[0,0,0,F(1,8)],[F(3,4),0,-3,F(3,4)]])
d = M.det()
print("det:", d, type(d))
print("det == F(9,512):", d == F(9,512))
print("F(str(det)) == F(9,512):", F(str(d)) == F(9,512))
