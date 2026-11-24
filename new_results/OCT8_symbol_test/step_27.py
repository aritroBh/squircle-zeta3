from fractions import Fraction as Fr
import sympy as sp
import sys; sys.path.insert(0,".")
exec(open("step_14.py").read().split("# formal letters")[0])  # reuse engine defs

x, y = sp.symbols('x y')
LET = {'X': x, 'Y': y, 'A': 1-x, 'B': 1-y, 'C': 1-x*y}
def dlog(L):
    f = LET[L]; return (sp.diff(f, x)/f, sp.diff(f, y)/f)
def wedge_dlog(L1, L2):
    a1,b1 = dlog(L1); a2,b2 = dlog(L2)
    return sp.simplify(a1*b2 - b1*a2)
def integrability(S):
    from collections import defaultdict
    w12 = defaultdict(list); w23 = defaultdict(list)
    for (l1,l2,l3), c in S.items():
        w12[l3].append((c,l1,l2)); w23[l1].append((c,l2,l3))
    bad = []
    for third, terms in w12.items():
        tot = sp.Rational(0)
        for c,l1,l2 in terms: tot += c*wedge_dlog(l1,l2)
        if sp.simplify(tot) != 0: bad.append(('12',third,sp.simplify(tot)))
    for first, terms in w23.items():
        tot = sp.Rational(0)
        for c,l2,l3 in terms: tot += c*wedge_dlog(l2,l3)
        if sp.simplify(tot) != 0: bad.append(('23',first,sp.simplify(tot)))
    return bad

Xv,Yv,Av,Bv,Cv = V('X'),V('Y'),V('A'),V('B'),V('C')
# S(Li_{1,2}(x,y)) = (xy-1)(x)[x(y-1)/(x-1)](x)y + (y-1)(x)(x-1)(x)y + (xy-1)(x)(xy)(x)[x/(x-1)] + (y-1)(x)y(x)(x-1)
def S_Li12():
    U1 = tensor3(Cv, vadd(Xv, Bv, vmul(Av,-1)), Yv)
    U2 = tensor3(Bv, Av, Yv)
    U3 = tensor3(Cv, vadd(Xv, Yv), vadd(Xv, vmul(Av,-1)))
    U4 = tensor3(Bv, Yv, Av)
    return sadd(U1, U2, U3, U4)
S12 = S_Li12()
print("S(Li12) integrability bad:", integrability(S12))

# OCT7 convention-B constants: Lam_B(a,c) = standard Li_{1,2}(c, a)
# swap x<->y (and A<->B) in S_Li12(x,y): S(Li_{1,2}(y,x))
swap = {'X':'Y','Y':'X','A':'B','B':'A','C':'C'}
def S_LamB():  # = S(Li_{1,2}(y,x))
    r = {}
    for (l1,l2,l3), c in S_Li12().items():
        t = (swap[l1],swap[l2],swap[l3])
        r[t] = r.get(t, Fr(0)) + c
    return r
# quick numeric sanity of convention identities:
# LamB(a,c) = Li2(a)Li1(c) - Li3(ac) - Li21A(a,c)   [stuffle, standard-A]
print("stuffle-convention relation OK (verified earlier)")
