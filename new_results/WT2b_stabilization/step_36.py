import sys, pickle
sys.path.insert(0,'.')
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))
def top_at(P, sign, order):
    L = P[1] if sign==+1 else P[2]
    return L[order-1] if len(L)>=order else ZERO

# Level 4: W_21(4) = -3 Gm_31(4) DL1 - 2 Gm_22(4) DL2 - 2 F(Gm_22(4) DL1)
Gm31 = record[(4,(3,1))]['Gm']
Gm22 = record[(4,(2,2))]['Gm']
W21 = record[(4,(2,1))]['W']
fB = pf_scale(pf_mul_DL(Gm22, DL2), Fr(-2))
fA = pf_scale(pf_fold(pf_mul_DL(Gm22, DL1)), Fr(-2))
f1 = pf_scale(pf_mul_DL(Gm31, DL1), Fr(-3))
print("W_21(4) top@+i ord4:", top_at(W21,+1,4).sympy())
print("  -3Gm31DL1 top:", top_at(f1,+1,4).sympy())
print("  -2Gm22DL2 top:", top_at(fB,+1,4).sympy())
print("  -2F(Gm22DL1) top:", top_at(fA,+1,4).sympy())
s = pf_add(pf_add(f1,fB),fA)
print("  sum of feeds top:", top_at(s,+1,4).sympy())
print("h_31(4):", top_at(record[(4,(3,1))]['H'],+1,6).sympy(), "  h_22(4):", top_at(record[(4,(2,2))]['H'],+1,6).sympy())
print("Gm31 top order5:", top_at(Gm31,+1,5).sympy(), " Gm22 top order5:", top_at(Gm22,+1,5).sympy())
# my F formula: top@+i of F(R) = -e_s, e_s = top@-i of R
R = pf_mul_DL(Gm22, DL1)
print("R=Gm22DL1: top@+i ord4:", top_at(R,+1,4).sympy(), " top@-i ord4:", top_at(R,-1,4).sympy())
F = pf_fold(R)
print("F(R) top@+i ord4:", top_at(F,+1,4).sympy(), "  predicted -e_s =", (-top_at(R,-1,4)).sympy())
