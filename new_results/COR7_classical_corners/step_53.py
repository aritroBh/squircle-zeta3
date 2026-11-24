import sys, json
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from fractions import Fraction as Fr
from cor6_sym import *
from cor6_table import vecQ

Q2 = W4Quotient(2)
def li4_sym(x, co=Fr(1)):
    vx, vo = vecQ(x), vecQ(1 - x)
    return {(lneg(vo), vx, vx, vx): co}

# COR6 T4 class at c=1/3
ST4 = {}
for x, co in ((Fr(-8), Fr(-11, 81)), (Fr(-3), Fr(-2)), (Fr(2, 3), Fr(20, 3))):
    ST4 = sadd(ST4, li4_sym(x, co))
cls_T4 = Q2.coords(ST4)
print("T4 class coords:", cls_T4)

# my numerical F4 Li4 part at c=1/3
Snum = {}
for x, co in ((Fr(-2), Fr(-22, 3)), (Fr(1, 3), Fr(-17, 2)), (Fr(2, 3), Fr(20, 3)),
              (Fr(1, 2), Fr(110, 9)), (Fr(1, 9), Fr(25, 48))):
    Snum = sadd(Snum, li4_sym(x, co))
cls_num = Q2.coords(Snum)
print("numerical F4 Li4-part class coords:", cls_num)
print("EQUAL:", cls_T4 == cls_num)

# also the COR6 omega symbol at c=1/3 itself
SOm3 = {}
for e, w, val in __import__('cor6_pipe').W4TERMS:
    SOm3 = sadd(SOm3, sscale(val * Fr(1, 3)**e, sym_H(w, Fr(1, 3), vecQ)))
print("class(Omega(1/3)) coords:", Q2.coords(SOm3))
print("Omega == T4:", Q2.coords(SOm3) == cls_T4)
