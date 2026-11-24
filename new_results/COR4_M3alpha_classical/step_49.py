import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from cor4_sym import sunit_vec_sq2
from fractions import Fraction as Fr
h = Fr(1,2)
print(sunit_vec_sq2((Fr(1), -h)), "expect (-1/2,-1) for 1-1/sqrt2")
print(sunit_vec_sq2((Fr(1), h)), "expect (-1/2,1) for 1+1/sqrt2")
print(sunit_vec_sq2((Fr(-1), Fr(1))), "expect (0,-1) for b=sqrt2-1")
print(sunit_vec_sq2((Fr(3), Fr(-2))), "expect (3/2,-1) for 1+b^2 = 4-2sqrt2? no: (3,-2)=3-2sqrt2=b^2 -> (0,-2)")
print(sunit_vec_sq2((Fr(2), Fr(-1))), "expect (1/2,-1) for 2-sqrt2")
print(sunit_vec_sq2((Fr(4), Fr(-2))), "expect (3/2,-1) for 4-2sqrt2 = 1+b^2")
