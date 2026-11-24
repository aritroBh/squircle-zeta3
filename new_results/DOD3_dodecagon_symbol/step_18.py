import sys, os
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

# ===== Validation of symbol formulas at FORMAL-letter level (exact), mirroring OCT8 step_14 =====
def V_(l): return {l: Fr(1)}
x, y = sp.symbols('x y')
LET = {'X': x, 'Y': y, 'A': 1-x, 'B': 1-y, 'C': 1-x*y}
def dlog(L):
    f = LET[L]
    return (sp.diff(f,x)/f, sp.diff(f,y)/f)
def wedge_dlog(L1,L2):
    a1,b1 = dlog(L1); a2,b2 = dlog(L2)
    return sp.simplify(a1*b2 - b1*a2)
def integrability(S):
    from collections import defaultdict
    w12 = defaultdict(list); w23 = defaultdict(list)
    for (l1,l2,l3),c in S.items():
        w12[l3].append((c,l1,l2)); w23[l1].append((c,l2,l3))
    bad12={}; bad23={}
    for third,terms in w12.items():
        tot = sp.Rational(0)
        for c,l1,l2 in terms: tot += c*wedge_dlog(l1,l2)
        if sp.simplify(tot)!=0: bad12[third]=sp.simplify(tot)
    for first,terms in w23.items():
        tot = sp.Rational(0)
        for c,l2,l3 in terms: tot += c*wedge_dlog(l2,l3)
        if sp.simplify(tot)!=0: bad23[first]=sp.simplify(tot)
    return bad12,bad23

Xv,Yv,Av,Bv,Cv = V_('X'),V_('Y'),V_('A'),V_('B'),V_('C')
def tensor3f(u,v,w,coeff=Fr(1)):
    S={}
    for a,ca in u.items():
        for b,cb in v.items():
            for c,cc in w.items():
                S[(a,b,c)]=S.get((a,b,c),Fr(0))+coeff*ca*cb*cc
    return {k:val for k,val in S.items() if val!=0}
def saddf(*Ss):
    r={}
    for S in Ss:
        for k,c in S.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def smulf(S,a):
    a=Fr(a); return {k:c*a for k,c in S.items() if c*a!=0}

# S(Li_{2,1}) formal
T1 = tensor3f(Cv, vadd(Xv,Yv), Bv)
T2 = tensor3f(Cv, vadd(Xv,Bv,vmul(Av,-1)), Xv)
T3 = tensor3f(Bv, Av, Xv)
S21 = saddf(T1,T2,T3)
b12,b23 = integrability(S21)
print("S(Li_{2,1}) #terms:", len(S21))
print("integrability wedge12 residual:", b12)
print("integrability wedge23 residual:", b23)

# stuffle: Li2(x)Li1(y) = Li21(x,y)+Li12(y,x)+Li3(xy)
def S_Li12f():
    U1 = tensor3f(Cv, vadd(Xv,Bv,vmul(Av,-1)), Yv)
    U2 = tensor3f(Bv, Av, Yv)
    U3 = tensor3f(Cv, vadd(Xv,Yv), vadd(Xv,vmul(Av,-1)))
    U4 = tensor3f(Bv, Yv, Av)
    return saddf(U1,U2,U3,U4)
swap = {'X':'Y','Y':'X','A':'B','B':'A','C':'C'}
S12_yx = {}
for (l1,l2,l3),c in S_Li12f().items():
    S12_yx[(swap[l1],swap[l2],swap[l3])] = S12_yx.get((swap[l1],swap[l2],swap[l3]),Fr(0))+c
S_Li3xy = tensor3f(vmul(Cv,-1), vadd(Xv,Yv), vadd(Xv,Yv))
S_stuffle = {}
for t in [('A','X','B'),('A','B','X'),('B','A','X')]: S_stuffle[t]=Fr(1)
res = saddf(S21, S12_yx, S_Li3xy, smulf(S_stuffle,-1))
print("stuffle residual (must be 0):", res)

# Landen for Li3: Li3(x)+Li3(1-x)+Li3(x/(x-1)) = products
def S_Li3f(omz, z): return tensor3f(vmul(omz,-1), z, z)
L1 = S_Li3f(Av, Xv)
L2 = S_Li3f(Xv, Av)
L3 = S_Li3f(vmul(Av,-1), vadd(Xv,vmul(Av,-1)))
prod = {}
for t in [('A','A','X'),('A','X','A'),('X','A','A')]: prod[t]=Fr(-1)
prod[('A','A','A')]=Fr(1)
res2 = saddf(L1,L2,L3, smulf(prod,-1))
print("Landen residual (must be 0):", res2)

# ===== Quotient over F =====
Q = Quotient()
print("\nDecomp rank over F (expect 19):", Q.rank_D, " dim Q (expect 8):", 27-Q.rank_D)

# ===== numeric spot check: S_Li3 at an F-point matches vecfun computation
z = P(2,-1,2)   # x_6 = (2-sqrt3)/2
Sz = S_Li3(z)
print("S(Li3(x_6)) terms:", len(Sz), " e.g.", dict(list(Sz.items())[:4]))
# direct: -(1-z) x z x z with 1-z = sqrt3/2
zv = vecfun(z); ozv = vecfun(psub(ONE,z))
print("  vec(x_6) =", zv, " vec(1-x_6) =", ozv)
man = tensor3(vmul(ozv,-1), zv, zv)
print("  matches manual:", Sz==man)
