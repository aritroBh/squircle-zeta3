import sympy as sp, pickle, os, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
x1 = -1/phi; x4 = 1/phi**2; kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
lf = sp.Symbol('lf')
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
def p1norm(e):
    return e.replace(lambda f: isinstance(f, sp.Function) and f.func == sp.polylog and f.args[0]==1,
                     lambda f: -sp.log(1-f.args[1]))
def cert(F, f, name):
    d = sp.diff(F, t) - f
    d = p1norm(sp.expand(d))
    # simplify with log(1+kap t) forms: expand logs
    d2 = sp.simplify(d)
    ok = (d2 == 0)
    print(f"cert {name}: {'OK' if ok else 'FAIL: '+str(d2)}")
    return ok
# F5 = phi[lf*w1 - Li2(x4*(1+t/phi))]  for w4/u1
F5 = phi*(lf*w1 - Li2(x4*(1+t/phi)))
# need d/dt of w4 and the log(1-x4(1+kap1 t)) = w4 - lf identity; sympy may not simplify directly; do manual check:
cert(F5, w4/u1, "F5 (w4/u1)")
