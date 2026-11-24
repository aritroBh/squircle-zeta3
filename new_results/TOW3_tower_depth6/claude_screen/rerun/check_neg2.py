import sys; sys.path.insert(0,'.')
import sympy as sp
import tow3lib as T3, prove_V56 as PV, prove_adj56 as PA
p = PA.p
A4 = T3.A[4]; T3.A[4] = A4 + 1
print('A_4 5->6:  V(6,2) numerator==0 ?', PV.identity(6, 2)); T3.A[4] = A4
A3 = T3.A[3]; T3.A[3] = A3 + 1
print('A_3 2->3:  V(5,1) numerator==0 ?', PV.identity(5, 1)); T3.A[3] = A3
def adj_resid(dl, s, pv):
    lhs = sp.expand(PA.eR(dl, s, PA.d, p) + PA.eR(dl, s+1, PA.d, p))
    return sp.simplify((lhs - sum(PA.feed_terms(dl, s))).subs(p, pv))
A3 = T3.A[3]; T3.A[3] = A3 + 1
print('A_3 2->3: adj (5,1) residual p=3,5:', [adj_resid(5, 1, pv) for pv in (3, 5)]); T3.A[3] = A3
A4 = T3.A[4]; T3.A[4] = A4 + 1
print('A_4 5->6: adj (6,2) residual p=3,5:', [adj_resid(6, 2, pv) for pv in (3, 5)]); T3.A[4] = A4
print('restored adj (6,2) residual p=3,5:', [adj_resid(6, 2, pv) for pv in (3, 5)])
