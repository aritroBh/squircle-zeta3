
"""Generalized OCT6 mixed engine for an edge with rational non-local ratios.
Term values: A0*t, A1*(1-t), B1*(1+k1 t), B2*(1+k2 t);  x_j = -k_j, beta_j = k_j/(1+k_j).
Named constants (symbols): l2,l3,l5 (logs), z2,z3, Li2/Li3 at
  x1,x2,b1,b2,-1/2,-2/3,-3/2,-2,3/5,2/3  and Lam1..Lam4 (depth-2).
"""
import pickle, sympy as sp
from functools import lru_cache

t = sp.symbols('t')
k1 = sp.Rational(1,3); k2 = sp.Rational(1,9)
x1 = -k1; x2 = -k2
b1 = k1/(1+k1); b2 = k2/(1+k2)     # 1/4, 1/10
ell1 = 2*sp.Symbol('l2') - sp.Symbol('l3')          # log(4/3)
ell2 = sp.Symbol('l2') + sp.Symbol('l5') - 2*sp.Symbol('l3')  # log(10/9)
l2, l3, l5, z2, z3 = sp.symbols('l2 l3 l5 z2 z3')
# single polylog atoms
L2x1, L3x1 = sp.symbols('L2x1 L3x1')   # Li2(-1/3), Li3(-1/3)
L2x2, L3x2 = sp.symbols('L2x2 L3x2')   # Li2(-1/9), Li3(-1/9)
L2b1, L3b1 = sp.symbols('L2b1 L3b1')   # Li2(1/4),  Li3(1/4)
L2b2, L3b2 = sp.symbols('L2b2 L3b2')   # Li2(1/10), Li3(1/10)
L2mh        = sp.symbols('L2mh')       # Li2(-1/2)
L2m2t       = sp.symbols('L2m2t')      # Li2(-2/3)
L2m3h, L3m3h= sp.symbols('L2m3h L3m3h')# Li2(-3/2), Li3(-3/2)
L2m2, L3m2  = sp.symbols('L2m2 L3m2')  # Li2(-2),   Li3(-2)
L2_3f, L3_3f= sp.symbols('L2_3f L3_3f')# Li2(3/5),  Li3(3/5)
L2_2t, L3_2t= sp.symbols('L2_2t L3_2t')# Li2(2/3),  Li3(2/3)
Lam1, Lam2, Lam3, Lam4 = sp.symbols('Lam1 Lam2 Lam3 Lam4')
NAMED = [l2,l3,l5,z2,z3,L2x1,L3x1,L2x2,L3x2,L2b1,L3b1,L2b2,L3b2,L2mh,L2m2t,
         L2m3h,L3m3h,L2m2,L3m2,L2_3f,L3_3f,L2_2t,L3_2t,Lam1,Lam2,Lam3,Lam4]

# ---- local table ----
import os
Jloc = pickle.load(open("research/08_theory/OCT5_octagon_exact/Jloc.pkl","rb"))
def JLOC(n,a,b):
    e = Jloc[(n,a,b)]
    return sp.expand(e.subs(sp.pi**2, 6*z2).subs(sp.zeta(3), z3))

# ---- K1d single-ratio table (symbolic in atoms at x) ----
K1d="research/08_theory/K1d_engine_fixed"
Jtab = pickle.load(open(K1d+'/Jtab_t.pkl','rb'))
L1s,L2s,L3s,AAs,z2s,xs = sp.symbols('L1 L2 L3 AA z2 x')
# Koelbig: A(x) = S_{1,2}(x) = 2Li3(x) + Li3(x/(x-1)) - Li2(x)*log(1-x) - log(1-x)^3/6
AA1 = 2*L3x1 + L3b1 - L2x1*ell1 - ell1**3/6   # x1/(x1-1) = 1/4 = b1, log(1-x1)=ell1
AA2 = 2*L3x2 + L3b2 - L2x2*ell2 - ell2**3/6   # x2/(x2-1) = 1/10 = b2
sub1 = {L1s: -ell1, L2s: L2x1, L3s: L3x1, AAs: AA1, z2s: z2, xs: x1}
sub2 = {L1s: -ell2, L2s: L2x2, L3s: L3x2, AAs: AA2, z2s: z2, xs: x2}
def K1(n,a,b,c,side):
    return sp.expand(Jtab[(n,a,b,c)].subs(sub1 if side==1 else sub2, simultaneous=True))

# ---- D constants (verified by quadrature, step_02) ----
D1 = -Lam1/k1
D2 = -Lam2/k2
D3 = sp.expand((Lam3/b1 - ell2*L2b1/b1)/(1+k1))
D4 = sp.expand((Lam4/b2 - ell1*L2b2/b2)/(1+k2))
rho = k2/(k1-k2); R1 = 1+k1
log_23 = l2 - l3     # log(2/3) = log((k1-k2)/k1), as symbols
D5 = sp.expand((log_23*ell1 + L2mh - L2m2t)/k1)
D6 = sp.expand((ell1*ell2 - k1*D5)/k2)
u0 = 1 - k1/k2   # -2
log_52 = l5 - l2     # log(5/2) = log(1-u0/R1)
D7 = sp.expand((ell1**2*log_52 - 2*ell1*L2m3h - 2*L3m3h + ell1**3/3 + 2*L3m2)/k2)
u0b = 1 - k2/k1  # 2/3
R2 = 1+k2
log_25 = l2 - l5     # log(2/5) = log(1-u0b/R2)
D8 = sp.expand((ell2**2*log_25 - 2*ell2*L2_3f - 2*L3_3f + ell2**3/3 + 2*L3_2t)/k1)
D9  = sp.expand((ell1**2*ell2 - k2*D7)/(2*k1))
D10 = sp.expand((ell1*ell2**2 - k1*D8)/(2*k2))

@lru_cache(None)
def MixedW(n):
    bnd = ell1*ell2
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
    return sp.expand((+ell1**2*ell2 - 2*k1*E1 - k2*E2)/(n+1))
@lru_cache(None)
def M4(n):
    F1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,0,2,2) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D8)
    F2 = sp.expand(sum((-1/k2)**j * MixedW(n-j) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D10)
    return sp.expand((ell1*ell2**2 - k1*F1 - 2*k2*F2)/(n+1))
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
