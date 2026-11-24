import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
x1 = -1/phi; x4 = 1/phi**2; kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
lf = sp.log(phi)
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)

def log_lin_key(arg):
    """arg = c*(1+e*t) (positive on [0,1]); return (c, e) with e simplified."""
    arg = sp.expand(arg)
    c = sp.simplify(arg.subs(t,0))
    rem = sp.simplify(arg - c)
    e = sp.simplify(rem/t) if rem != 0 else sp.Integer(0)
    return c, sp.expand(e)

def cert2(F, f, name):
    D = sp.diff(F, t) - f
    D = D.replace(lambda g: isinstance(g, sp.Function) and g.func == sp.polylog and g.args[0]==1,
                  lambda g: -sp.log(1-g.args[1]))
    D = sp.expand(D)
    # collect log terms
    logterms = {}   # key e -> coeff sum
    constlogs = []  # (coeff, c)
    rest = 0
    for term in sp.Add.make_args(D):
        dd = term.as_powers_dict()
        logs = [b for b in dd if isinstance(b, sp.log)]
        if not logs:
            rest += term
            continue
        if len(logs) != 1 or dd[logs[0]] != 1:
            print("non-simple log term:", term); return False
        arg = logs[0].args[0]
        coef = sp.simplify(term/logs[0])
        c, e = log_lin_key(arg)
        if e == 0:
            constlogs.append((coef, c))
        else:
            logterms.setdefault(e, []).append(coef)
    ok = True
    r = sp.simplify(rest)
    if r != 0:
        print(f"  {name}: rational part nonzero:", r); ok = False
    for e, cl in logterms.items():
        s = sp.simplify(sum(cl))
        if s != 0:
            print(f"  {name}: log(1+({e})t) coeff nonzero:", s); ok = False
    if constlogs:
        prod = sp.Integer(1)
        for coef, c in constlogs:
            prod = prod * c**coef
        prod = sp.simplify(prod)
        if prod != 1:
            print(f"  {name}: constant-log product != 1:", prod); ok = False
    print(f"cert {name}: {'OK' if ok else 'FAIL'}")
    return ok

F5 = phi*(lf*w1 - Li2(x4*(1+t/phi)))
cert2(F5, w4/u1, "F5 (D5: w4/u1)")
F6 = -phi**2*(2*lf*w4 - Li2((1/phi)*(1+kap4*t)))
cert2(F6, w1/u4, "F6 (D6: w1/u4)")
s7 = 1+kap4*t
G7 = 4*lf**2*sp.log(s7) - 4*lf*Li2(s7/phi) + sp.log(s7/phi)*sp.log(1-s7/phi)**2 \
     + 2*sp.log(1-s7/phi)*Li2(1-s7/phi) - 2*Li3(1-s7/phi)
F7 = -phi**2*G7
cert2(F7, w1**2/u4, "F7 (D7: w1^2/u4)")
s8 = 1+t/phi
H8 = lf**2*sp.log(s8) - 2*lf*Li2(x4*s8) + sp.log(x4*s8)*sp.log(1-x4*s8)**2 \
     + 2*sp.log(1-x4*s8)*Li2(1-x4*s8) - 2*Li3(1-x4*s8)
F8 = phi*H8
cert2(F8, w4**2/u1, "F8 (D8: w4^2/u1)")
