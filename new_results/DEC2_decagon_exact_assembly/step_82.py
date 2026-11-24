import sympy as sp, pickle, os, sys, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
kap1 = 1/phi; kap4 = -1/phi**2
import dec_engine as E
Jtab = pickle.load(open('Jtab_dec2.pkl','rb'))
D = pickle.load(open('D_all_proved.pkl','rb'))
def K1v(n,a,b,c,side):
    return sp.expand(Jtab[(n,a,b,c)].subs(E.sub1 if side==1 else E.sub4, simultaneous=True))
def canon5(e):
    e = sp.expand(e)
    P = sp.Poly(e, *E.NAMED)
    out = sp.Integer(0)
    for monom, coeff in P.terms():
        c = sp.together(coeff)
        cc = c.subs(r5, -r5)
        p = sp.radsimp(sp.expand((c+cc)/2)); q = sp.radsimp(sp.expand((c-cc)/(2*r5)))
        mm = sp.Integer(1)
        for a_, e_ in zip(E.NAMED, monom): mm *= a_**e_
        out += (p + q*r5)*mm
    return sp.expand(out)
from functools import lru_cache
@lru_cache(None)
def MixedW(n):
    bnd = -E.lf**2
    A1 = sp.expand(sum((-1/kap1)**j * K1v(n-j,0,0,1,4) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*D['D5'])
    A2 = sp.expand(sum((-1/kap4)**j * K1v(n-j,0,0,1,1) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*D['D6'])
    return sp.expand((bnd - kap1*A1 - kap4*A2)/(n+1))
@lru_cache(None)
def M1(n):
    B1 = sp.expand(sum((-1/kap1)**j * K1v(n-j,1,0,1,4) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*D['D1'])
    B2 = sp.expand(sum((-1/kap4)**j * K1v(n-j,1,0,1,1) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*D['D2'])
    return sp.expand(-(MixedW(n) + kap1*B1 + kap4*B2)/(n+1))
@lru_cache(None)
def M2(n):
    r1 = (-1/kap1)**(n+1) - 1
    r4 = (-1/kap4)**(n+1) - 1
    C1 = sp.expand(sum((-1/kap1)**j * K1v(n-j,0,1,1,4) for j in range(n+1))/kap1 + r1*D['D3'])
    C2 = sp.expand(sum((-1/kap4)**j * K1v(n-j,0,1,1,1) for j in range(n+1))/kap4 + r4*D['D4'])
    return sp.expand((-sum(MixedW(j) for j in range(n+1)) - kap1*C1 - kap4*C2)/(n+1))
@lru_cache(None)
def M3(n):
    E1 = sp.expand(sum((-1/kap1)**j * MixedW(n-j) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*D['D9'])
    E2 = sp.expand(sum((-1/kap4)**j * K1v(n-j,0,0,2,1) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*D['D7'])
    return sp.expand((-E.lf**3 - 2*kap1*E1 - kap4*E2)/(n+1))
@lru_cache(None)
def M4(n):
    F1 = sp.expand(sum((-1/kap1)**j * K1v(n-j,0,0,2,4) for j in range(n+1))/kap1 + (-1/kap1)**(n+1)*D['D8'])
    F2 = sp.expand(sum((-1/kap4)**j * MixedW(n-j) for j in range(n+1))/kap4 + (-1/kap4)**(n+1)*D['D10'])
    return sp.expand((E.lf**3 - kap1*F1 - 2*kap4*F2)/(n+1))
def MIXm(n,a,b,c1,c4):
    if a==1: return M1(n)
    if b==1: return M2(n)
    if (c1,c4)==(2,1): return M3(n)
    if (c1,c4)==(1,2): return M4(n)
    return MixedW(n)
t0=time.time()
bad = 0
MIXvals = {}
for typ in [(0,0,1,1),(1,0,1,1),(0,1,1,1),(0,0,2,1),(0,0,1,2)]:
    for n in range(4):
        a,b,c1,c4 = typ
        mine = MIXm(n,a,b,c1,c4)
        eng = E.MIX(n,a,b,c1,c4)
        d = canon5(mine - eng)
        MIXvals[(n,a,b,c1,c4)] = mine
        if d != 0:
            bad += 1
            print("MISMATCH", (n,a,b,c1,c4), str(d)[:200])
print("mixed integrals computed: 20, symbolic mismatches vs DEC1 engine:", bad, " %.0fs"%(time.time()-t0))
pickle.dump(MIXvals, open('MIX20_proved.pkl','wb'))
