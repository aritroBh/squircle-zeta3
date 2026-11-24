import mpmath as mp
mp.mp.dps = 200
phi = (1+mp.sqrt(5))/2
lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi**2
L31 = mp.polylog(3, 1/phi); L31m = mp.polylog(3, -1/phi); L34 = mp.polylog(3, 1/phi**2)
def pslq_test(name, val, basis, bnames, tol=1e-150):
    import sympy as sp
    ints = [mp.mpf(v) for v in basis] + [val]
    # simple LLL-based PSLQ via sympy? use mp.pslq? mpmath has pslq
    rel = mp.pslq(ints, tol=tol, maxcoeff=10**12, maxsteps=10000)
    print(name, "->", rel)
    if rel:
        # check
        r = sum(c*v for c, v in zip(rel, ints))
        print("   residual:", mp.nstr(abs(r),3), " relation:", dict(zip(bnames+['TARGET'], rel)))
print("T1: Li3(1/phi) in span{z3, pi2*lf, lf^3}?")
pslq_test("Li3(1/phi)", L31, [z3, pi2*lf, lf**3], ['z3','pi2*lf','lf^3'])
print("T2: Li3(-1/phi)?")
pslq_test("Li3(-1/phi)", L31m, [z3, pi2*lf, lf**3], ['z3','pi2*lf','lf^3'])
print("T3: Li3(1/phi^2)?")
pslq_test("Li3(1/phi^2)", L34, [z3, pi2*lf, lf^3 if False else lf**3], ['z3','pi2*lf','lf^3'])
print("T4: Li4(-1/phi) in span{pi4, z3*lf, pi2*lf^2, lf^4}?")
pslq_test("Li4(-1/phi)", mp.polylog(4,-1/phi), [mp.pi**4, z3*lf, pi2*lf**2, lf**4], ['pi4','z3lf','pi2lf2','lf4'])
print("T4b: Li4(1/phi^2)?")
pslq_test("Li4(1/phi^2)", mp.polylog(4,1/phi**2), [mp.pi**4, z3*lf, pi2*lf**2, lf**4], ['pi4','z3lf','pi2lf2','lf4'])
print("T5: Li5(-1/phi) in span{z5, pi2*z3, z3*lf^2, pi2*lf^3, pi4*lf, lf^5}?")
pslq_test("Li5(-1/phi)", mp.polylog(5,-1/phi), [z5, pi2*z3, z3*lf**2, pi2*lf**3, mp.pi**4*lf, lf**5], ['z5','pi2z3','z3lf2','pi2lf3','pi4lf','lf5'])
print("T5b: Li5(1/phi^2)?")
pslq_test("Li5(1/phi^2)", mp.polylog(5,1/phi**2), [z5, pi2*z3, z3*lf**2, pi2*lf**3, mp.pi**4*lf, lf**5], ['z5','pi2z3','z3lf2','pi2lf3','pi4lf','lf5'])
