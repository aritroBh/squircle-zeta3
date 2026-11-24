import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
x1 = -1/phi; x4 = 1/phi**2; kap1 = 1/phi; kap4 = -1/phi**2
lf = sp.log(phi); l2 = sp.log(2)
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
def unit2phi(c):
    c = sp.simplify(c)
    for b in range(-8,9):
        v2 = sp.simplify(c/phi**b)
        q = sp.simplify((v2 - v2.subs(r5,-r5))/(2*r5))
        if q == 0:
            p = sp.simplify((v2 + v2.subs(r5,-r5))/2)
            try: pr = sp.Rational(p)
            except TypeError: continue
            fa = sp.factorint(int(pr.p)); fb = sp.factorint(int(pr.q))
            if set(fa)|set(fb) <= {2}:
                a = fa.get(2,0)-fb.get(2,0)
                assert sp.simplify(c - 2**a*phi**b) == 0
                return (a,b)
    return None
def canon_logs(e):
    letters = {}
    def repl(f):
        arg = sp.expand(f.args[0])
        c = sp.simplify(arg.subs(t,0))
        rem = sp.expand(arg - c)
        ev = sp.simplify(rem/(c*t)) if rem != 0 else sp.Integer(0)
        ev = sp.expand(ev)
        ab = unit2phi(c)
        if ab is None: raise ValueError(f"bad log arg: {arg} (c={c})")
        base = ab[0]*l2 + ab[1]*lf
        if ev == 0:
            return base
        if ev not in letters:
            letters[ev] = sp.log(1+ev*t)
        return base + letters[ev]
    e2 = e.replace(lambda f: isinstance(f, sp.log), repl)
    return e2, letters
def cert3(F, f, name):
    D = sp.diff(F, t) - f
    D = D.replace(lambda g: isinstance(g, sp.Function) and g.func == sp.polylog and g.args[0]==1,
                  lambda g: -sp.log(1-g.args[1]))
    D, letters = canon_logs(D)
    D = sp.expand(D)
    lfs, l2s = sp.Symbol('LFs'), sp.Symbol('L2S')
    D2 = D.subs({lf: lfs, l2: l2s})
    gens = list(letters.values()) + [lfs, l2s]
    P = sp.Poly(sp.expand(D2), *gens)
    bad = 0
    for monom, coeff in P.terms():
        c = sp.simplify(coeff)
        if c != 0:
            bad += 1
            print(f"  {name}: NONZERO coeff of {monom}: {str(c)[:150]}")
    print(f"cert {name}: {'OK (formal identity, all coeffs 0)' if bad==0 else 'FAIL'}")
    return bad==0
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
F5 = phi*(lf*w1 - Li2(x4*(1+t/phi)))
c5 = cert3(F5, w4/u1, "F5 (D5: w4/u1)")
F6 = -phi**2*(2*lf*w4 - Li2((1/phi)*(1+kap4*t)))
c6 = cert3(F6, w1/u4, "F6 (D6: w1/u4)")
s7 = 1+kap4*t
G7 = 4*lf**2*sp.log(s7) - 4*lf*Li2(s7/phi) + sp.log(s7/phi)*sp.log(1-s7/phi)**2 \
     + 2*sp.log(1-s7/phi)*Li2(1-s7/phi) - 2*Li3(1-s7/phi)
F7 = -phi**2*G7
c7 = cert3(F7, w1**2/u4, "F7 (D7: w1^2/u4)")
s8 = 1+t/phi
H8 = lf**2*sp.log(s8) - 2*lf*Li2(x4*s8) + sp.log(x4*s8)*sp.log(1-x4*s8)**2 \
     + 2*sp.log(1-x4*s8)*Li2(1-x4*s8) - 2*Li3(1-x4*s8)
F8 = phi*H8
c8 = cert3(F8, w4**2/u1, "F8 (D8: w4^2/u1)")
print("all pass:", c5 and c6 and c7 and c8)
