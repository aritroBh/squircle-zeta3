import sympy as sp, pickle, os, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
x1 = -1/phi; x4 = 1/phi**2; kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
lf = sp.log(phi)   # use actual log(phi) now for certificates
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
def norm_log(e):
    e = sp.expand_log(e, force=True)
    e = e.replace(lambda f: isinstance(f, sp.Function) and f.func == sp.polylog and f.args[0]==1,
                  lambda f: -sp.log(1-f.args[1]))
    e = sp.expand_log(e, force=True)
    return sp.simplify(sp.expand(e))
def cert(F, f, name):
    d = norm_log(sp.diff(F, t)) - norm_log(f)
    d = sp.simplify(d)
    ok = (d == 0)
    print(f"cert {name}: {'OK' if ok else 'FAIL: '+str(d)[:300]}")
    return ok
# F5 for w4/u1
F5 = phi*(lf*w1 - Li2(x4*(1+t/phi)))
cert(F5, w4/u1, "F5 (D5: w4/u1)")
# F6 for w1/u4 : F6 = -phi^2 [2 lf w4 - Li2((1/phi)(1+kap4 t))]
F6 = -phi**2*(2*lf*w4 - Li2((1/phi)*(1+kap4*t)))
cert(F6, w1/u4, "F6 (D6: w1/u4)")
# F7 for w1^2/u4 : F7 = -phi^2 G(s), s = 1+kap4 t
s7 = 1+kap4*t
G7 = 4*lf**2*sp.log(s7) - 4*lf*Li2(s7/phi) + sp.log(s7/phi)*sp.log(1-s7/phi)**2 \
     + 2*sp.log(1-s7/phi)*Li2(1-s7/phi) - 2*Li3(1-s7/phi)
F7 = -phi**2*G7
cert(F7, w1**2/u4, "F7 (D7: w1^2/u4)")
# F8 for w4^2/u1 : F8 = phi H(s), s = 1+t/phi
s8 = 1+t/phi
H8 = lf**2*sp.log(s8) - 2*lf*Li2(x4*s8) + sp.log(x4*s8)*sp.log(1-x4*s8)**2 \
     + 2*sp.log(1-x4*s8)*Li2(1-x4*s8) - 2*Li3(1-x4*s8)
F8 = phi*H8
cert(F8, w4**2/u1, "F8 (D8: w4^2/u1)")
