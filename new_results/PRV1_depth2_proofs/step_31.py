
part1 = r'''

# ======================================================================
# Part 1: classical dilog/trilog function identities on 0<t<1 (exact)
# ======================================================================
# Formal differentiation engine.  Letters: R=log t, A=log(1-t), B=log(1+t),
# C=log 2, Z=zeta(2), Z3=zeta(3).  Polylog symbols carry a derivative map
# obtained from d/dz Li_n(z) = Li_{n-1}(z)/z with real branches on (0,1).
t = sp.symbols('t', positive=True)
R, A, B, C, Z, Z3 = sp.symbols('R A B C Z Z3')
Dt, Dm, Db, Dh, Dc, Dmc, D1, Du, Dw = sp.symbols('Dt Dm Db Dh Dc Dmc D1 Du Dw')
Dt2 = sp.Symbol('Dt2')
Tt, Tm, Tmc, Tu, T1, Tb, Tw, Tc = sp.symbols('Tt Tm Tmc Tu T1 Tb Tw Tc')
Tt2 = sp.Symbol('Tt2')
q, lc, U = B-C, A-B, C+R-B        # log((1+t)/2), log c, log u ; c=(1-t)/(1+t), u=2t/(1+t)
DERIV = {
    R: 1/t, A: -1/(1-t), B: 1/(1+t), C: 0, Z: 0, Z3: 0,
    Dt: -A/t, Dm: -B/t, D1: R/(1-t), Db: q/(1-t), Dh: -(A-C)/(1+t),
    Dc: 2*U/(1-t**2), Dmc: -2*q/(1-t**2), Dw: B/(t*(1+t)),
    Dt2: -2*(A+B)/t,
    Tt: Dt/t, Tm: Dm/t, Tmc: -2*Dmc/(1-t**2), Tu: Du/(t*(1+t)),
    T1: -D1/(1-t), Tb: -Db/(1-t), Tw: Dw/(t*(1+t)), Tc: -2*Dc/(1-t**2),
    Tt2: 2*Dt2/t,
}
diff = lambda e: sum(sp.diff(e, v)*dv for v, dv in DERIV.items())

def part1():
    print("=" * 74)
    print("PART 1: classical dilog/trilog function identities (derivative = 0 exactly)")
    ids = {
        "E   Euler: Li2(t)+Li2(1-t)=Z-log t log(1-t)": Dt + D1 + R*A - Z,
        "L1  Landen at -t: Li2(-t)+Li2(w)=-B^2/2":     Dm + Dw + B**2/2,
        "L2  Landen at -c: Li2(-c)+Li2(b)=-q^2/2":     Dmc + Db + q**2/2,
        "DUP2 dilog duplication: Li2(t^2)=2Li2(t)+2Li2(-t)": Dt2 - 2*Dt - 2*Dm,
        "OC  odd complement: Li2(t)-Li2(-t)+Li2(c)-Li2(-c)=3Z/2-log t log c":
            Dt - Dm + Dc - Dmc + R*lc - 3*Z/2,
        "DUP3 trilog duplication: Li3(t^2)=4Li3(t)+4Li3(-t)": Tt2 - 4*Tt - 4*Tm,
    }
    for name, e in ids.items():
        r = sp.factor(diff(e))
        print("   d/dt[", name, "] =", r)
        assert r == 0
    # DUP3's derivative uses DUP2 as a function identity: check that reduction too
    r = sp.factor(diff(Tt2 - 4*Tt - 4*Tm).subs(Dt2, 2*Dt+2*Dm))
    print("   d/dt[DUP3] after DUP2 =", r)
    assert r == 0
    print("  exact constant pins:")
    print("   E  at t->0+: Li2(0)+Li2(1)=Z ; log t log(1-t)->0   => 0  (exact)")
    print("   L1 at t=0:   0+0+0 = 0  (exact)")
    print("   L2 at t=0:   Li2(-1)+Li2(1/2)+log^2(1/2)/2")
    print("              = -Z/2 + (Z/2 - C^2/2) + C^2/2 = 0")
    print("      [Li2(-1)=-Z/2 from DUP2 at t=1; Li2(1/2)=Z/2-C^2/2 from E at t=1/2]")
    print("   DUP2 at t=0: 0  (exact);   OC at t->1-: Z+Z/2-3Z/2=0, log t log c->0 (exact)")
    print("   DUP3 at t=0: 0  (exact);   Li3(-1) = -3*Z3/4 follows from DUP3 at t=1.")
    print("PART 1 OK")
'''
with open('verify_full.py','a') as f: f.write(part1)
print("part1 appended")
