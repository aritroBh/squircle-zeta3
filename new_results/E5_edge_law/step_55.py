import mpmath as mp, json
mp.mp.dps = 100
P5 = mp.mpf(json.load(open('steps/step_10_num.json'))['P5'])
z2, z3, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(5)
V = P5/mp.sqrt(5)
basis = [mp.mpf(1), z2, z3, z2**2, z2*z3, z5]
bn = ['1','z2','z3','z2^2','z2z3','z5']
for tol in ['1e-45', '1e-50', '1e-55']:
    rel = mp.pslq(mp.matrix([V]+basis), tol=mp.mpf(tol), maxcoeff=10**8, maxsteps=10000)
    print("tol", tol, "->", rel)
    if rel is not None and rel[0] != 0:
        from fractions import Fraction
        q = [Fraction(int(rel[i+1]), int(-rel[0])) for i in range(6)]
        print("   V =", {bn[i]: str(q[i]) for i in range(6) if q[i] != 0})
        resid = V + sum(mp.mpf(rel[i+1])*basis[i]/rel[0] for i in range(6))
        print("   residual:", mp.nstr(abs(resid), 3))
        break
