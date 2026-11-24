import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
t = sp.symbols('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
k1 = 1/phi; k4 = -1/phi**2
x1 = -k1; x4 = -k4
l2, lf, z2, z3 = sp.symbols('l2 lf z2 z3')
L2_1, L3_1, AA_1 = sp.symbols('L2_1 L3_1 AA_1')   # Li2(x1), Li3(x1), A(x1)
L2_4, L3_4, AA_4 = sp.symbols('L2_4 L3_4 AA_4')   # at x4
L2f, L3f = sp.symbols('L2f L3f')                  # Li2(1/phi), Li3(1/phi)
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
# D constants (decagon, derived & quad-verified)
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
# save engine state
eng = {'NAMED':NAMED,'D':{i:sp.expand(v) for i,v in enumerate([None,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10]) if v is not None},'k1':k1,'k4':k4,'x1':x1,'x4':x4,'phi':phi}
pickle.dump(eng, open('engine_state.pkl','wb'))
# quick check: evaluate a few mixed integrals numerically
mp.mp.dps = 110
v = pickle.load(open('atoms830.pkl','rb'))
phim = v['phi']; lfm=v['lf']; l2m=v['l2']; z2m=v['z2']; z3m=v['z3']
num = {l2:l2m, lf:lfm, z2:z2m, z3:z3m,
       L2_1: mp.polylog(2,-1/phim), L3_1: v['L3_x1'], AA_1: v['A_x1'],
       L2_4: mp.polylog(2,1/phim**2), L3_4: v['L3_x4'], AA_4: v['A_x4'],
       L2f: mp.polylog(2,1/phim), L3f: v['L3_xf'], Lam: v['Lam'], swap: v['swap']}
def ev(e): return mp.mpf(sp.N(e.subs(num, simultaneous=True), 110))
kap1m=1/phim; kap4m=-1/phim**2
w1m=lambda t: mp.log(1+kap1m*t); w4m=lambda t: mp.log(1+kap4m*t)
tests = []
for n in range(4):
    tests += [(M1(n), lambda t,n=n: t**n*mp.log(t)*w1m(t)*w4m(t), f'M1({n})'),
              (M2(n), lambda t,n=n: t**n*mp.log(1-t)*w1m(t)*w4m(t), f'M2({n})'),
              (M3(n), lambda t,n=n: t**n*w1m(t)**2*w4m(t), f'M3({n})'),
              (M4(n), lambda t,n=n: t**n*w1m(t)*w4m(t)**2, f'M4({n})'),
              (MixedW(n), lambda t,n=n: t**n*w1m(t)*w4m(t), f'MixedW({n})')]
worst=mp.mpf(0)
for e,f,nm in tests:
    d = abs(ev(e)-mp.quad(f,[0,1])); worst=max(worst,d)
    if d>mp.mpf('1e-95'): print("BAD", nm, mp.nstr(d,4))
print("worst |exact - quad| over 20 mixed integrals (110 dps):", mp.nstr(worst,4))
with open('progress.log','a') as f:
    f.write(f"step: decagon mixed engine built (D1..D10 closed forms, recursions). 20 mixed integrals vs quad 110dps: worst {mp.nstr(worst,4)}\n")