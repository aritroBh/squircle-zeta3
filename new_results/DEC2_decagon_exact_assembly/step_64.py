import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
x1 = -1/phi; x4 = 1/phi**2; kap1 = 1/phi; kap4 = -1/phi**2
lf = sp.log(phi); l2 = sp.log(2)
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)

def unit2phi(c):
    """c in Q(sqrt5)^+ : return (a,b) with c = 2^a * phi^b, or None."""
    c = sp.simplify(c)
    for b in range(-8,9):
        v = sp.expand(c/phi**b)
        v2 = sp.simplify(v)
        # check rational
        q = sp.simplify((v2 - v2.subs(r5,-r5))/(2*r5))
        if q == 0:
            p = sp.simplify((v2 + v2.subs(r5,-r5))/2)
            try:
                pr = sp.Rational(p)
            except TypeError:
                continue
            from sympy import factorint
            fa = sp.factorint(int(pr.p)); fb = sp.factorint(int(pr.q))
            if set(fa)|set(fb) <= {2}:
                a = fa.get(2,0)-fb.get(2,0)
                assert sp.simplify(c - 2**a*phi**b) == 0
                return (a,b)
    return None

def canon_logs(e):
    """rewrite every log(arg) with arg = c*(1+e t) as a*l2 + b*lf + log(1+e t); return (expr, letter set)."""
    letters = {}
    def repl(f):
        arg = sp.expand(f.args[0])
        c = sp.simplify(arg.subs(t,0))
        rem = sp.simplify(arg - c)
        ev = sp.expand(rem/t) if rem != 0 else sp.Integer(0)
        ab = unit2phi(c)
        if ev == 0:
            if ab is None: raise ValueError(f"constant log not 2^a phi^b: {c}")
            return ab[0]*l2 + ab[1]*lf
        if ab is None: raise ValueError(f"bad linear log arg: {arg}")
        key = ev
        if key not in letters:
            letters[key] = sp.log(1+ev*t)
        return ab[0]*l2 + ab[1]*lf + letters[key]
    e2 = e.replace(lambda f: isinstance(f, sp.log), repl)
    return e2, letters

def cert3(F, f, name):
    D = sp.diff(F, t) - f
    D = D.replace(lambda g: isinstance(g, sp.Function) and g.func == sp.polylog and g.args[0]==1,
                  lambda g: -sp.log(1-g.args[1]))
    D, letters = canon_logs(D)
    D = sp.expand(D)
    # letter inventory: log-letters + polylogs + {lf, l2} as symbols
    lfs, l2s = sp.Symbol('LFs'), sp.Symbol('L2S')
    D2 = D.subs({lf: lfs, l2: l2s})
    gens = list(letters.values()) + [lfs, l2s]
    polys = []
    for g in letters.values():
        pass
    # collect coefficients wrt products of letter functions
    P = sp.Poly(sp.expand(D2), *gens) if gens else sp.Poly(sp.expand(D2))
    bad = 0
    for monom, coeff in P.terms():
        c = sp.simplify(coeff)
        if c != 0:
            bad += 1
            print(f"  {name}: NONZERO coeff of {monom}: {str(c)[:200]}")
    print(f"cert {name}: {'OK (all coefficients 0 -> formal identity)' if bad==0 else 'FAIL'}")
    return bad==0

u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
F5 = phi*(lf*w1 - Li2(x4*(1+t/phi)))
cert3(F5, w4/u1, "F5 (D5: w4/u1)")
F6 = -phi**2*(2*lf*w4 - Li2((1/phi)*(1+kap4*t)))
cert3(F6, w1/u4, "F6 (D6: w1/u4)")
s7 = 1+kap4*t
G7 = 4*lf**2*sp.log(s7) - 4*lf*Li2(s7/phi) + sp.log(s7/phi)*sp.log(1-s7/phi)**2 \
     + 2*sp.log(1-s7/phi)*Li2(1-s7/phi) - 2*Li3(1-s7/phi)
F7 = -phi**2*G7
cert3(F7, w1**2/u4, "F7 (D7: w1^2/u4)")
s8 = 1+t/phi
H8 = lf**2*sp.log(s8) - 2*lf*Li2(x4*s8) + sp.log(x4*s8)*sp.log(1-x4*s8)**2 \
     + 2*sp.log(1-x4*s8)*Li2(1-x4*s8) - 2*Li3(1-x4*s8)
F8 = phi*H8
cert3(F8, w4**2/u1, "F8 (D8: w4^2/u1)")
