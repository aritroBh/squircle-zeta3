engine_src = '''
import pickle, sympy as sp
K1d="research/08_theory/K1d_engine_fixed"
t = sp.symbols('t'); r2 = sp.sqrt(2)
k1 = r2-1; k2 = r2/sp.Integer(2)-1
x1 = -k1; x2 = -k2
l2, lu, z2, z3 = sp.symbols('l2 lu z2 z3')
L2_1, L3_1, AA_1 = sp.symbols('L2_1 L3_1 AA_1')
L2_2, L3_2, AA_2 = sp.symbols('L2_2 L3_2 AA_2')
L2k, L3k, L2d, L3d = sp.symbols('L2k L3k L2d L3d')
Lam, swap = sp.symbols('Lam swap')
L1_1 = -l2/2; L1_2 = l2/2
NAMED = [l2,lu,z2,z3,L2_1,L3_1,AA_1,L2_2,L3_2,AA_2,L2k,L3k,L2d,L3d,Lam,swap]
import os
O5 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'oct5')
Jloc = pickle.load(open(os.path.join(O5,'Jloc.pkl'),'rb'))
def JLOC(n,a,b):
    e = Jloc[(n,a,b)]
    return sp.expand(e.subs(sp.pi**2, 6*z2).subs(sp.zeta(3), z3))
Jtab = pickle.load(open(K1d+'/Jtab_t.pkl','rb'))
L1s,L2s,L3s,AAs,z2s,xs = sp.symbols('L1 L2 L3 AA z2 x')
sub1 = {L1s:L1_1, L2s:L2_1, L3s:L3_1, AAs:AA_1, z2s:z2, xs:x1}
sub2 = {L1s:L1_2, L2s:L2_2, L3s:L3_2, AAs:AA_2, z2s:z2, xs:x2}
def K1(n,a,b,c,side):
    return sp.expand(Jtab[(n,a,b,c)].subs(sub1 if side==1 else sub2, simultaneous=True))
D1 = -Lam/k1
D2 = -swap/k2
D3 = sp.simplify(((l2/2)*L2_2/x2 + swap/x2)/r2)
D4 = sp.simplify((Lam/x1 - (l2/2)*L2_1/x1)/(r2/2))
D5 = sp.expand(((lu-l2/2)*(l2/2) + L2k - L2d)/k1)
D6 = sp.expand((-lu*l2/2 - L2k + L2d)/k2)
D7 = sp.expand((-(l2**2/4)*lu + l2*L2d - 2*L3d + 2*L3k)/k2)
D8 = sp.expand(((l2**2/4)*(l2/2-lu) - l2*L2k - 2*L3k + 2*L3d)/k1)
D9 = sp.expand(-l2**3/(16*k1) - (k2/(2*k1))*D7)
D10= sp.expand(l2**3/(16*k2) - (k1/(2*k2))*D8)
from functools import lru_cache
@lru_cache(None)
def MixedW(n):
    bnd = -l2**2/4
    A1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,0,1,2) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D5)
    A2 = sp.expand(sum((-1/k2)**j * K1(n-j,0,0,1,1) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D6)
    return sp.expand((bnd - k1*A1 - k2*A2)/(n+1))
@lru_cache(None)
def M1(n):
    B1 = sp.expand(sum((-1/k1)**j * K1(n-j,1,0,1,2) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D1)
    B2 = sp.expand(sum((-1/k2)**j * K1(n-j,1,0,1,1) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D2)
    return sp.expand(-(MixedW(n) + k1*B1 + k2*B2)/(n+1))
@lru_cache(None)
def M2(n):
    r1 = (-1/k1)**(n+1) - 1
    r2v = (-1/k2)**(n+1) - 1
    C1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,1,1,2) for j in range(n+1))/k1 + r1*D3)
    C2 = sp.expand(sum((-1/k2)**j * K1(n-j,0,1,1,1) for j in range(n+1))/k2 + r2v*D4)
    return sp.expand((-sum(MixedW(j) for j in range(n+1)) - k1*C1 - k2*C2)/(n+1))
@lru_cache(None)
def M3(n):
    E1 = sp.expand(sum((-1/k1)**j * MixedW(n-j) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D9)
    E2 = sp.expand(sum((-1/k2)**j * K1(n-j,0,0,2,1) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D7)
    return sp.expand((-l2**3/8 - 2*k1*E1 - k2*E2)/(n+1))
@lru_cache(None)
def M4(n):
    F1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,0,2,2) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D8)
    F2 = sp.expand(sum((-1/k2)**j * MixedW(n-j) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D10)
    return sp.expand((l2**3/8 - k1*F1 - 2*k2*F2)/(n+1))
def MIX(n,a,b,c1,c2):
    if a==1: return M1(n)
    if b==1: return M2(n)
    if (c1,c2)==(2,1): return M3(n)
    if (c1,c2)==(1,2): return M4(n)
    return MixedW(n)
def INT(n,a,b,c1,c2):
    if c1==0 and c2==0: return JLOC(n,a,b)
    if c2==0: return K1(n,a,b,c1,1)
    if c1==0: return K1(n,a,b,c2,2)
    return MIX(n,a,b,c1,c2)
'''
open('mixed_engine.py','w').write(engine_src)
print("mixed_engine.py written", len(engine_src), "chars")
