
import sys
sys.path.insert(0,'.')
from fastred import *
from fractions import Fraction as Fr
tau1 = QSI(-1,0,-1,0)/QSI(2); tau2 = QSI(-1,0,1,0)/QSI(2)

def mk_G(glist):
    # PF with +i principal part glist [g1..gM] and conjugate -i part
    P = pf_zero()
    P[1] = list(glist)
    P[2] = [e.conj_i() for e in glist]
    return P

ok1 = ok2 = ok3 = True
for dd in range(2, 14):
    M = 2*dd-1
    for trial in range(3):
        gv = [QSI(Fr(2*j+1+trial,j+2), Fr(j,j+3), Fr(j-1+trial,j+4), Fr(j,j+5)) for j in range(1, M+1)]
        G = mk_G(gv)
        prod = pf_mul_DL(G, DL1)
        tc = prod[1][2*dd-3] if len(prod[1])>=2*dd-2 else ZERO
        ok1 &= (tc == gv[M-1]*tau1)
        prod2 = pf_mul_DL(G, DL2)
        tc2 = prod2[1][2*dd-3] if len(prod2[1])>=2*dd-2 else ZERO
        ok2 &= (tc2 == gv[M-1]*tau2)
    # lower terms only:
    gv0 = [QSI(Fr(j+1,j+2),0,Fr(j,j+3),0) for j in range(1, M)] 
    G0 = mk_G(gv0)
    prod0 = pf_mul_DL(G0, DL1)
    tc0 = prod0[1][2*dd-3] if len(prod0[1])>=2*dd-2 else ZERO
    ok3 &= tc0.is_zero()
print("(B1) [top](G*DL1) = tau1_0*[top of G], d=2..13, exact fastred:", ok1)
print("(B1') same for DL2 with tau2_0:", ok2)
print("(B1'') lower terms of G never feed the top of the product:", ok3)

# B2: Int factor: HigherPart top h at order 2d -> Gm top at order 2d-1 is -h/(2d-1)
ok4 = True
for dd in range(2,14):
    h = QSI(Fr(3,7),Fr(1,5),Fr(2,9),Fr(1,4))
    H = pf_zero(); H[1]=[ZERO]*(2*dd-2); H[1].append(h)  # orders 1..2d-1? need h at order 2d
    H[1] = [ZERO]*(2*dd-1) + [h]   # h at order 2d
    H[2] = [e.conj_i() for e in H[1]]
    G = pf_integrate_higher(H)
    top = G[1][2*dd-2] if len(G[1])>=2*dd-1 else ZERO
    ok4 &= (top == -h/QSI(2*dd-1))
print("(B2) Int: top at order 2d maps to top at order 2d-1 with factor -1/(2d-1), d=2..13:", ok4)

# B3: F-branch factor -1 in fastred arithmetic
ok5 = True
for m in range(1, 30):
    c = QSI(Fr(2,5),Fr(1,3),Fr(3,7),Fr(1,2))
    R = pf_zero(); R[2]=[ZERO]*(m-1)+[c]   # pole at -i order m
    R[1]=[e.conj_i() for e in R[2]]
    FR = pf_fold(R)
    top = FR[1][m-1] if len(FR[1])>=m else ZERO
    ok5 &= (top == -c)
print("(B3) pf_fold: top at -i order m -> top at +i order m with factor -1, m=1..29:", ok5)
