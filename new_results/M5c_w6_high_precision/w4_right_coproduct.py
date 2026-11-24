"""(1,3) coaction check. Right-hand periods are taken modulo pi^2.

Primitives below extend the supplied two-log integrator to arbitrary letter
endpoints, including finite parts. Branch changes proportional to pi^2 are
discarded ONLY in the right de Rham factor, never in w4_coproduct.py.
"""
from w4_coproduct import *


def modpi(p):
    return P({k:c for k,c in p.d.items() if 'pi2' not in k})


def cube(p):
    return p*p*p


def hprimitive(c, w):
    """Primitive of log^2|w|/(w-c), modulo pi^2, unit tangent in w."""
    assert c != ZERO
    if w == ZERO:
        return P()
    lc = logp(c)
    if w == c:
        return -cube(lc)-PE.Z3*2
    s = w/c
    return modpi(lc*lc*logp(s-ONE) + lc*PE.P1(s)*2 + PE.P2(s))


def sqprimitive(x, z, e):
    if x == z:
        return cube(logp(e-x))*F(1,3)
    return hprimitive(z-x, e-x)


def singleprimitive(x, z, e):
    if x == z:
        l = logp(e-x)
        return l*l*F(1,2)
    lc = logp(z-x)
    if e == x:
        return P()
    if e == z:
        return -lc*lc
    s = (e-x)/(z-x)
    return modpi(lc*logp(s-ONE)+PE.P1(s))


def ratioprimitive(x, y, z, e):
    """Primitive of log^2|(v-x)/(v-y)|/(v-z), FP in v at v=e."""
    assert x != y
    if z == x:
        if e == x:
            return -cube(logp(x-y))*F(1,3)
        if e == y:
            return P()
        w = (e-x)/(e-y)
        return cube(logp(w))*F(1,3)-hprimitive(ONE,w)
    if z == y:
        if e == y:
            return -cube(logp(x-y))*F(1,3)
        if e == x:
            return P()
        return -hprimitive(ONE,(e-x)/(e-y))
    c = (x-z)/(y-z)
    if e == x:
        return P()
    if e == y:
        # H_c(w)-H_1(w) tends to -log^3|c|/3 as w tends to infinity.
        return -cube(logp(c))*F(1,3)
    w = (e-x)/(e-y)
    ans = hprimitive(c,w)-hprimitive(ONE,w)
    if e == z:
        # log|w-c| = log|v-z| + log|w'(z)| + o(1).
        lc = logp(c)
        ans += lc*lc*logp((x-y)/((z-y)*(z-y)))
    return ans


def pairprimitive(x, y, z, e):
    if x == y:
        return sqprimitive(x,z,e)
    return (sqprimitive(x,z,e)+sqprimitive(y,z,e)-ratioprimitive(x,y,z,e))*F(1,2)


@cache
def primitive(z, mixed, e):
    (d1,c1), (d2,c2) = LT_LOG, MT_LOG if mixed else LT_LOG
    out = P()
    for x,cx in d1.items():
        for y,cy in d2.items():
            out += pairprimitive(x,y,z,e)*cx*cy
    for x,cx in d1.items():
        out += singleprimitive(x,z,e)*cx*c2
    for y,cy in d2.items():
        out += singleprimitive(y,z,e)*cy*c1
    return modpi(out+logp(e-z)*c1*c2)


def interval(z, mixed, lo, hi):
    return primitive(z,mixed,hi)-primitive(z,mixed,lo)


def eval_log(lf,e):
    d,c = lf
    return c+sum((logp(e-x)*cx for x,cx in d.items()),P())


def left_tensor(ell, p):
    return [p*ell.d.get((name,),ZERO) for name in ('lam','L')]


def lhs_right_coproduct():
    def component(z,mixed):
        if not mixed:
            forms = [(LT_LOG,False,3)]
            pz = cube(eval_log(LT_LOG,z))
        else:
            forms = [(LT_LOG,True,2),(MT_LOG,False,1)]
            lz = eval_log(LT_LOG,z)
            pz = lz*lz*eval_log(MT_LOG,z)
        out = left_tensor(logp(A-z)-logp(B-z),pz)
        for lf, flag, multiplicity in forms:
            d,_ = lf
            for x,cx in d.items():
                left_a = logp(A-x)-logp(ONE-x)
                left_b = logp(B-x)-logp(ONE-x)
                out = vadd(out,vscale(left_tensor(left_a,interval(z,flag,x,A)),multiplicity*cx),
                           vscale(left_tensor(left_b,interval(z,flag,x,B)),-multiplicity*cx))
            out = vadd(out,vscale(left_tensor(eval_log(lf,ONE),interval(z,flag,B,A)),multiplicity))
        return [modpi(p) for p in out]
    return vadd(vscale(component(ZERO,False),(8-7*R)/48),
                vscale(component(-A,False),Q(F(-1,3))),
                vscale(component(ZERO,True),R/16))


def rhs_right_coproduct():
    import math
    out = [P(),P()]
    for c,atoms in conjecture_terms():
        if ('pi2',) in atoms:
            continue  # pi^2 is zero in the right factor and cannot have left weight 1
        for i,a in enumerate(atoms):
            if a[0] == 'log':
                da = tensor(P.const(1),a[1])
            elif a[0] == 'li':
                n,z = a[1:]
                right = P.const(1)
                for _ in range(n-1):
                    right *= logp(z)
                da = left_tensor(-logp(ONE-z),right*F(1,math.factorial(n-1)))
            else:
                assert a == ('z3',)
                continue
            for j,b in enumerate(atoms):
                if i != j:
                    da = [v*atom_p(b) for v in da]
            out = vadd(out,vscale(da,c))
    return out


def main():
    ban_float()
    import sympy as sp
    v,x,y,z = sp.symbols('v x y z')
    w,c = (v-x)/(v-y),(x-z)/(y-z)
    assert sp.factor(sp.diff(w,v)*(1/(w-c)-1/(w-1))-1/(v-z)) == 0
    assert sp.factor(-sp.diff(w,v)/(w-1)-1/(v-y)) == 0
    assert sp.factor(sp.diff(w,v)*(1/w-1/(w-1))-1/(v-x)) == 0
    W,C = sp.symbols('W C')
    assert sp.expand(C*C*(W-C)+C*(W-C)**2+(W-C)**3/3-(W**3-C**3)/3) == 0
    print('Three exact Mobius differential checks and infinity-limit algebra: PASS',flush=True)
    # Independent exact comparison with the existing two-log integration code.
    for z,flag in [(ZERO,False),(-A,False),(A,False),(-B,False),(ZERO,True),(-A,True),(B,False),(B,True)]:
        delta = PE.reduce_li2_all(interval(z,flag,B,A)-modpi(pair(z,flag)))
        for a in list(delta.atoms('Li3')):
            if a < -ONE:
                l = logp(-a)
                val = PE.LI3(ONE/a)-PE.PI2*l*F(1,6)-cube(l)*F(1,6)
                delta = delta.subs_atom(('Li3',a),val)
        if modpi(delta).d:
            delta = li3_reduce(delta)
        assert not modpi(delta).d, (z,flag,delta)
    print('Generalized primitive versus supplied two-log integrator: PASS',flush=True)
    lhs,rhs = lhs_right_coproduct(),rhs_right_coproduct()
    for i,name in enumerate(['log2','log a']):
        print('Reducing left coefficient of',name,flush=True)
        red = modpi(li3_reduce(PE.reduce_li2_all(lhs[i]-rhs[i])))
        print('Delta_(1,3)(W4-RHS), left coefficient of',name,'=',red,flush=True)
        assert not red.d
    print('PROVED-MODULO-REDUCED-COACTION-KERNEL: Delta_(1,3)(W4-RHS) = 0')


if __name__ == '__main__':
    main()
