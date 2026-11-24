import sys, pickle
sys.path.insert(0,'.')
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))

# At level k=3: decompose W_21 = R_21 - 3 Gm_31 DL1 - 2 Gm_22 DL2 - 2 F(Gm_22 DL1)
W21 = record[(3,(2,1))]['W']
R21_init = record[(3,(2,1))]  # W at initial level includes feeds from (2,2)!
# W_21 at k=3 = R_21 + feeds from Gm_22(3) (same level!) and Gm_31(3)=0 (deg 4 > 3)
Gm22 = record[(3,(2,2))]['Gm']
print("Gm_22(3) +i PP:", [str(e.sympy()) for e in Gm22[1]], " const:", Gm22[3].sympy())

def top_at(P, sign, order):
    L = P[1] if sign==+1 else P[2]
    return L[order-1] if len(L)>=order else ZERO

# feeds:
B = pf_scale(pf_mul_DL(Gm22, DL2), Fr(-2))          # -2 Gm_22 DL2
A = pf_scale(pf_fold(pf_mul_DL(Gm22, DL1)), Fr(-2))  # -2 F(Gm_22 DL1)
print("B-feed top@+i (ord4):", top_at(B,+1,4).sympy())
print("F-feed top@+i (ord4):", top_at(A,+1,4).sympy())
print("W_21(3) top@+i (ord4):", top_at(W21,+1,4).sympy())
# so R_21 top = W21 - B - A
R21top = top_at(W21,+1,4) - top_at(B,+1,4) - top_at(A,+1,4)
print("=> R_21 (initial) top@+i:", R21top.sympy())
# what did my rho formula give? rho_21 = -3 sqrt2 /4
print("formula rho_21 = -3*sqrt2/4 =", (-3*S2/4))
# check Gm_22 top order and value
print("h_22(3) top (order 6? Gm top 5): Gm22 order5:", top_at(Gm22,+1,5).sympy())
H22 = record[(3,(2,2))]['H']
print("H_22(3) ipos:", [str(e.sympy()) for e in H22[1]])
# Wait: key (2,2) has degree 4, not present at level 3! Let me check what's in record
print("record has (3,(2,2))?", (3,(2,2)) in record)
