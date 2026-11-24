
import pickle, sympy as sp
t = sp.symbols('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
k1 = 1/phi; k4 = -1/phi**2
x1 = -k1; x4 = -k4
l2, lf, z2, z3 = sp.symbols('l2 lf z2 z3')
L2_1, L3_1, AA_1 = sp.symbols('L2_1 L3_1 AA_1')
L2_4, L3_4, AA_4 = sp.symbols('L2_4 L3_4 AA_4')
L2f, L3f = sp.symbols('L2f L3f')
Lam, swap = sp.symbols('Lam swap')
L1_1 = -lf; L1_4 = lf
NAMED = [l2,lf,z2,z3,L2_1,L3_1,AA_1,L2_4,L3_4,AA_4,L2f,L3f,Lam,swap]
O5='research/08_theory/OCT6_octagon_closed/oct5'
Jloc = pickle.load(open(O5+'/Jloc.pkl','rb'))
def JLOC(n,a,b):
    e = Jloc[(n,a,b)]
    return sp.expand(e.subs(sp.pi**2, 6*z2).subs(sp.zeta(3), z3))
Jtab = pickle.load(open('research/08_theory/K1d_engine_fixed/Jtab_t.pkl','rb'))
L1s,L2s,L3s,AAs,z2s,xs = sp.symbols('L1 L2 L3 AA z2 x')
sub1 = {L1s:L1_1, L2s:L2_1, L3s:L3_1, AAs:AA_1, z2s:z2, xs:x1}
sub4 = {L1s:L1_4, L2s:L2_4, L3s:L3_4, AAs:AA_4, z2s:z2, xs:x4}
def K1(n,a,b,c,side):
    return sp.expand(Jtab[(n,a,b,c)].subs(sub1 if side==1 else sub4, simultaneous=True))
D1 = sp.expand(Lam/x1)
D2 = sp.expand(swap/x4)
D3 = sp.expand(phi*(lf*L2_4 + swap))
D4 = sp.expand(phi**2*(lf*L2_1 - Lam))
D5 = sp.expand(phi*(lf**2 - L2f + L2_4))
D6 = sp.expand(phi**2*(L2_4 - L2f + 2*lf**2))
D7 = sp.expand(2*phi**2*(lf**3 - lf*L2f + L3f - L3_4))
D8 = sp.expand(phi*(2*L3f - lf**3 - 2*lf*L2_4 - 2*L3_4))
D9 = sp.expand(-phi*lf**3/2 + D7/(2*phi))
D10= sp.expand(-phi**2*lf**3/2 + phi*D8/2)
from functools import lru_cache
@lru_cache(None)
def MixedW(n):
    bnd = -lf**2
    A1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,0,1,4) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D5)
    A2 = sp.expand(sum((-1/k4)**j * K1(n-j,0,0,1,1) for j in range(n+1))/k4 + (-1/k4)**(n+1)*D6)
    return sp.expand((bnd - k1*A1 - k4*A2)/(n+1))
@lru_cache(None)
def M1(n):
    B1 = sp.expand(sum((-1/k1)**j * K1(n-j,1,0,1,4) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D1)
    B2 = sp.expand(sum((-1/k4)**j * K1(n-j,1,0,1,1) for j in range(n+1))/k4 + (-1/k4)**(n+1)*D2)
    return sp.expand(-(MixedW(n) + k1*B1 + k4*B2)/(n+1))
@lru_cache(None)
def M2(n):
    r1 = (-1/k1)**(n+1) - 1
    r4 = (-1/k4)**(n+1) - 1
    C1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,1,1,4) for j in range(n+1))/k1 + r1*D3)
    C2 = sp.expand(sum((-1/k4)**j * K1(n-j,0,1,1,1) for j in range(n+1))/k4 + r4*D4)
    return sp.expand((-sum(MixedW(j) for j in range(n+1)) - k1*C1 - k4*C2)/(n+1))
@lru_cache(None)
def M3(n):
    E1 = sp.expand(sum((-1/k1)**j * MixedW(n-j) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D9)
    E2 = sp.expand(sum((-1/k4)**j * K1(n-j,0,0,2,1) for j in range(n+1))/k4 + (-1/k4)**(n+1)*D7)
    return sp.expand((-lf**3 - 2*k1*E1 - k4*E2)/(n+1))
@lru_cache(None)
def M4(n):
    F1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,0,2,4) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D8)
    F2 = sp.expand(sum((-1/k4)**j * MixedW(n-j) for j in range(n+1))/k4 + (-1/k4)**(n+1)*D10)
    return sp.expand((lf**3 - k1*F1 - 2*k4*F2)/(n+1))
def MIX(n,a,b,c1,c4):
    if a==1: return M1(n)
    if b==1: return M2(n)
    if (c1,c4)==(2,1): return M3(n)
    if (c1,c4)==(1,2): return M4(n)
    return MixedW(n)
def INT(n,a,b,c1,c4):
    if c1==0 and c4==0: return JLOC(n,a,b)
    if c4==0: return K1(n,a,b,c1,1)
    if c1==0: return K1(n,a,b,c4,4)
    return MIX(n,a,b,c1,c4)
