import cor8r_sym as S
import sympy as sp
cc = S.cc
lv = -(1+cc)/(2*cc)
print("factor_Qc(-(1+c)/(2c)) =", S.factor_Qc(lv), "  expect (-1,-1,0,1) = [2c/(1+c)]^{-1}")
lv2 = (1+cc)/(2*cc)
print("factor_Qc((1+c)/(2c))  =", S.factor_Qc(lv2))
lv3 = sp.cancel((1-2)/(2-2/(1+cc)))
print("letter value (1-2)/(2-B) =", lv3)
print("factor_Qc of it        =", S.factor_Qc(lv3))
