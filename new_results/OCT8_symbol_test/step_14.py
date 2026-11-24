from fractions import Fraction as Fr
import sympy as sp

# ---------- symbol engine over formal letters ----------
# A 'vec' is dict letter->Fraction (element of V = F^x (x) Q, letters formal).
# A 'symbol' is dict (l1,l2,l3)->Fraction with letters expanded (pure tensors over letter basis).

def vadd(*vs):
    r = {}
    for v in vs:
        for k, c in v.items():
            r[k] = r.get(k, Fr(0)) + c
            if r[k] == 0: del r[k]
    return r
def vmul(v, a):
    a = Fr(a)
    return {k: c*a for k, c in v.items() if c*a != 0}
def V(letter): return {letter: Fr(1)}

def tensor3(u, v, w, coeff=Fr(1)):
    # u,v,w vecs -> symbol dict
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            for c, cc in w.items():
                S[(a,b,c)] = S.get((a,b,c), Fr(0)) + coeff*ca*cb*cc
    return {k: x for k, x in S.items() if x != 0}
def sadd(*Ss):
    r = {}
    for S in Ss:
        for k, c in S.items():
            r[k] = r.get(k, Fr(0)) + c
            if r[k] == 0: del r[k]
    return r
def smul(S, a):
    a = Fr(a)
    return {k: c*a for k, c in S.items() if c*a != 0}
def shuffle13(a, S2):
    # a letter, S2 = {(l1,l2):c}: a (sh) (l1 l2)
    r = {}
    for (l1,l2), c in S2.items():
        for t in [(a,l1,l2),(l1,a,l2),(l1,l2,a)]:
            r[t] = r.get(t, Fr(0)) + c
    return {k: x for k, x in r.items() if x != 0}

# formal letters as sympy functions of x,y for integrability
x, y = sp.symbols('x y')
LET = {'X': x, 'Y': y, 'A': 1-x, 'B': 1-y, 'C': 1-x*y}
def dlog(L):
    f = LET[L]
    return (sp.diff(f, x)/f, sp.diff(f, y)/f)  # (coeff of dx, coeff of dy)
def wedge_dlog(L1, L2):
    a1, b1 = dlog(L1); a2, b2 = dlog(L2)
    return sp.simplify(a1*b2 - b1*a2)  # coeff of dx^dy

def integrability(S):
    """return (bad12, bad23): dicts third-factor->2form (should all be 0)"""
    from collections import defaultdict
    w12 = defaultdict(list); w23 = defaultdict(list)
    for (l1,l2,l3), c in S.items():
        w12[l3].append((c, l1, l2))
        w23[l1].append((c, l2, l3))
    bad12 = {}
    for third, terms in w12.items():
        tot = sp.Rational(0)
        for c, l1, l2 in terms:
            tot += c*wedge_dlog(l1, l2)
        if sp.simplify(tot) != 0: bad12[third] = sp.simplify(tot)
    bad23 = {}
    for first, terms in w23.items():
        tot = sp.Rational(0)
        for c, l2, l3 in terms:
            tot += c*wedge_dlog(l2, l3)
        if sp.simplify(tot) != 0: bad23[first] = sp.simplify(tot)
    return bad12, bad23

# letter shorthands
Xv, Yv, Av, Bv, Cv = V('X'), V('Y'), V('A'), V('B'), V('C')

# S(Li_3(z)) = -(1-z) (x) z (x) z
def S_Li3(one_minus_z, z):
    return tensor3(vmul(one_minus_z, -1), z, z)

# S(Li_{2,1}(x,y)) = (xy-1)(x)(xy)(x)(y-1) + (xy-1)(x)[x(y-1)/(x-1)](x)x + (y-1)(x)(x-1)(x)x
def S_Li21():
    T1 = tensor3(Cv, vadd(Xv, Yv), Bv)
    T2 = tensor3(Cv, vadd(Xv, Bv, vmul(Av, -1)), Xv)
    T3 = tensor3(Bv, Av, Xv)
    return sadd(T1, T2, T3)

S21 = S_Li21()
print("S(Li_{2,1}) terms:", len(S21))
b12, b23 = integrability(S21)
print("Integrability wedge12 residual:", b12, " wedge23 residual:", b23)

# Stuffle: S(Li2 x * Li1 y) = (A(x)X) sh (B)
S_stuffle = shuffle13('A', {('X','B'): Fr(1)})  # wait: (A (x) X) sh B = A X B + A B X + B A X
S_stuffle = {}
for t in [('A','X','B'), ('A','B','X'), ('B','A','X')]:
    S_stuffle[t] = Fr(1)

# S(Li_{1,2}(y,x)): swap X<->Y, A<->B of S(Li_{1,2}(x,y))
def S_Li12():
    U1 = tensor3(Cv, vadd(Xv, Bv, vmul(Av,-1)), Yv)
    U2 = tensor3(Bv, Av, Yv)
    U3 = tensor3(Cv, vadd(Xv, Yv), vadd(Xv, vmul(Av,-1)))
    U4 = tensor3(Bv, Yv, Av)
    return sadd(U1, U2, U3, U4)
S12_yx = {}
swap = {'X':'Y','Y':'X','A':'B','B':'A','C':'C'}
for (l1,l2,l3), c in S_Li12().items():
    S12_yx[(swap[l1],swap[l2],swap[l3])] = S12_yx.get((swap[l1],swap[l2],swap[l3]),Fr(0))+c
S_Li3xy = tensor3(vmul(Cv,-1), vadd(Xv,Yv), vadd(Xv,Yv))
lhs = sadd(S21, S12_yx, S_Li3xy)
print("Stuffle check: S(Li21)+S(Li12(y,x))+S(Li3(xy)) - (A(x)X)sh B =", sadd(lhs, smul(S_stuffle,-1)))

# Landen for Li_3: S(Li3(x))+S(Li3(1-x))+S(Li3(x/(x-1))) =? S(products)
L1 = S_Li3(Av, Xv)                    # Li3(x): 1-z = 1-x = A, z = X
L2 = S_Li3(Xv, Av)                    # Li3(1-x)
L3 = S_Li3(vmul(Av,-1), vadd(Xv, vmul(Av,-1)))  # Li3(x/(x-1)): 1-z = 1/(1-x) ~ A^{-1}, z = x/(x-1) ~ X - A
landen_lhs = sadd(L1, L2, L3)
prod = {}
for t in [('A','A','X'),('A','X','A'),('X','A','A')]: prod[t] = Fr(-1)
prod[('A','A','A')] = Fr(1)
print("Landen check:", sadd(landen_lhs, smul(prod, -1)))
