import pickle, sympy as sp
from o6util import load

t = sp.symbols('t'); r2 = sp.sqrt(2)
k1 = r2-1; k2 = r2/sp.Integer(2)-1
x1 = -k1; x2 = -k2
# named constant symbols
l2, lu, z2, z3 = sp.symbols('l2 lu z2 z3')
L2_1, L3_1, AA_1 = sp.symbols('L2_1 L3_1 AA_1')
L2_2, L3_2, AA_2 = sp.symbols('L2_2 L3_2 AA_2')
L2k, L3k, L2d, L3d = sp.symbols('L2k L3k L2d L3d')
Lam, swap = sp.symbols('Lam swap')
L1_1 = -l2/2; L1_2 = l2/2   # -log(1-x1) = -l2/2 ; -log(1-x2) = +l2/2

# ---- local table
Jloc = load('Jloc.pkl')
pi, zeta3 = sp.symbols('pi zeta3')
def JLOC(n,a,b):
    return Jloc[(n,a,b)].xreplace({sp.pi**2: 6*z2, sp.pi: sp.sqrt(6*z2)}).xreplace({sp.zeta(3): z3})
def JLOC(n,a,b):
    e = Jloc[(n,a,b)]
    e = e.subs(sp.pi**2, 6*z2).subs(sp.zeta(3), z3)
    return sp.expand(e)

# ---- K1d single-ratio table
K1d="research/08_theory/K1d_engine_fixed"
Jtab = pickle.load(open(K1d+'/Jtab_t.pkl','rb'))
L1s,L2s,L3s,AAs,z2s,xs = sp.symbols('L1 L2 L3 AA z2 x')
sub1 = {L1s:L1_1, L2s:L2_1, L3s:L3_1, AAs:AA_1, z2s:z2, xs:x1}
sub2 = {L1s:L1_2, L2s:L2_2, L3s:L3_2, AAs:AA_2, z2s:z2, xs:x2}
def K1(n,a,b,c,side):
    return sp.expand(Jtab[(n,a,b,c)].subs(sub1 if side==1 else sub2, simultaneous=True))

# ---- D constants (exact closed forms, step_02 verified)
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
    # int t^n W1 W2 dt
    bnd = -l2**2/4
    A1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,0,1,2) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D5)
    A2 = sp.expand(sum((-1/k2)**j * K1(n-j,0,0,1,1) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D6)
    return sp.expand((bnd - k1*A1 - k2*A2)/(n+1))

@lru_cache(None)
def M1(n):  # int t^n X W1 W2
    B1 = sp.expand(sum((-1/k1)**j * K1(n-j,1,0,1,2) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D1)
    B2 = sp.expand(sum((-1/k2)**j * K1(n-j,1,0,1,1) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D2)
    return sp.expand(-(MixedW(n) + k1*B1 + k2*B2)/(n+1))

@lru_cache(None)
def M2(n):  # int t^n Y W1 W2
    r1 = (-1/k1)**(n+1) - 1
    r2v = (-1/k2)**(n+1) - 1
    C1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,1,1,2) for j in range(n+1))/k1 + r1*D3)
    C2 = sp.expand(sum((-1/k2)**j * K1(n-j,0,1,1,1) for j in range(n+1))/k2 + r2v*D4)
    return sp.expand((-sum(MixedW(j) for j in range(n+1)) - k1*C1 - k2*C2)/(n+1))

@lru_cache(None)
def M3(n):  # int t^n W1^2 W2
    E1 = sp.expand(sum((-1/k1)**j * MixedW(n-j) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D9)
    E2 = sp.expand(sum((-1/k2)**j * K1(n-j,0,0,2,1) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D7)
    return sp.expand((-l2**3/8 - 2*k1*E1 - k2*E2)/(n+1))

@lru_cache(None)
def M4(n):  # int t^n W1 W2^2
    F1 = sp.expand(sum((-1/k1)**j * K1(n-j,0,0,2,2) for j in range(n+1))/k1 + (-1/k1)**(n+1)*D8)
    F2 = sp.expand(sum((-1/k2)**j * MixedW(n-j) for j in range(n+1))/k2 + (-1/k2)**(n+1)*D10)
    return sp.expand((l2**3/8 - k1*F1 - 2*k2*F2)/(n+1))

def MIX(n,a,b,c1,c2):
    if (a,b,c1,c2)[0]==1: return M1(n)
    if b==1: return M2(n)
    if (c1,c2)==(2,1): return M3(n)
    if (c1,c2)==(1,2): return M4(n)
    return MixedW(n)

pickle.dump('ok', open('tmp_flag.pkl','wb'))
# stash the functions' source is in this step file; test numerically next block
import mpmath as mp
mp.mp.dps=110
r2m=mp.sqrt(2); k1m=r2m-1; k2m=r2m/2-1
num = {
 l2: mp.log(2), lu: mp.log(1+r2m), z2: mp.pi**2/6, z3: mp.zeta(3),
 L2_1: mp.polylog(2,-k1m), L3_1: mp.polylog(3,-k1m),
 L2_2: mp.polylog(2,-k2m), L3_2: mp.polylog(3,-k2m),
 L2k: mp.polylog(2,k1m), L3k: mp.polylog(3,k1m),
 L2d: mp.polylog(2,2-r2m), L3d: mp.polylog(3,2-r2m),
 Lam: mp.mpf(pickle.load(open('oct5/atoms600.pkl','rb'))['Li21_a']),
 swap: mp.mpf(pickle.load(open('oct5/atoms600.pkl','rb'))['Li21_b']),
}
def AA_of(xv,N=4000):
    s=mp.mpf(0); Hn=mp.mpf(0)
    for nn in range(1,N+1):
        Hn+=1/mp.mpf(nn); s+=xv**nn*Hn/nn**2
    return s
num[AA_1]=AA_of(-k1m); num[AA_2]=AA_of(-k2m)
def ev(e): return mp.mpf(sp.N(e.subs(num, simultaneous=True), 100))
print("numeric checks of mixed families vs quadrature (110 dps):")
tests = []
for n in range(4):
    tests += [(M1(n), lambda t,n=n: t**n*mp.log(t)*mp.log(1+k1m*t)*mp.log(1+k2m*t)),
              (M2(n), lambda t,n=n: t**n*mp.log(1-t)*mp.log(1+k1m*t)*mp.log(1+k2m*t)),
              (M3(n), lambda t,n=n: t**n*mp.log(1+k1m*t)**2*mp.log(1+k2m*t)),
              (M4(n), lambda t,n=n: t**n*mp.log(1+k1m*t)*mp.log(1+k2m*t)**2),
              (MixedW(n), lambda t,n=n: t**n*mp.log(1+k1m*t)*mp.log(1+k2m*t))]
worst=mp.mpf(0)
for e,f in tests:
    d = abs(ev(e)-mp.quad(f,[0,1])); worst=max(worst,d)
print("worst |exact - quad| over 20 mixed integrals:", mp.nstr(worst,3))
