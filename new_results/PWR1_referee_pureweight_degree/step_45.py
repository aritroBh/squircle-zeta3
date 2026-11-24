import mpmath as mp, json, sympy as sp
mp.mp.dps = 100
js = json.load(open('../D16_pure_weight/tau_series_2d.json'))
Z = {k: sp.Symbol(f'Z{k}') for k in range(2, 13)}
fit = [mp.mpf(s) for s in json.load(open('fit_w2_w12.json'))]
for k in [11, 12]:
    expr = 0
    for mono, c in js['logP'][str(k)].items():
        term = sp.Rational(c)
        for zi in mono.split():
            term *= Z[int(zi)]
        expr += term
    val = mp.mpf(str(expr.subs({Z[i]: sp.zeta(i) for i in Z}).evalf(50)))
    rel = abs(val-fit[k-2])/abs(val)
    print(f"w_{k}: exact={mp.nstr(val,25)} fit={mp.nstr(fit[k-2],25)} relerr={mp.nstr(rel,2)} ({mp.nstr(-mp.log10(rel),3)} digits)")

# PSLQ purity tests on fitted coefficients (weight-4 and weight-5)
zeta = mp.zeta
print("\nPSLQ purity: w4 should be a Q-combo of {z2^2, z4} ONLY (weight 4); test against bigger basis")
basis4 = [mp.mpf(1), mp.euler, zeta(2), zeta(3), zeta(2)**2, zeta(4)]
rel = mp.pslq([fit[2]] + basis4, tol=mp.mpf('1e-11'), maxcoeff=10000, maxsteps=10000)
print("  pslq(w4, 1, gamma, z2, z3, z2^2, z4):", rel)
# expected relation: w4 = -z2^2/4 + z4/32  -> 32*w4 + 8*z2^2 - z4 = 0
print("  check 32*w4 + 8*z2^2 - z4 =", mp.nstr(32*fit[2] + 8*zeta(2)**2 - zeta(4), 6))
print("\nPSLQ purity: w5 should be Q-combo of {z2 z3, z5} ONLY (weight 5)")
basis5 = [mp.mpf(1), mp.euler, zeta(2), zeta(3), zeta(2)**2, zeta(4), zeta(2)*zeta(3), zeta(5)]
rel = mp.pslq([fit[3]] + basis5, tol=mp.mpf('1e-9'), maxcoeff=10000, maxsteps=10000)
print("  pslq(w5, 1, gamma, z2, z3, z2^2, z4, z2z3, z5):", rel)
print("  check 16*w5 + 8*z2z3 + 3*z5 =", mp.nstr(16*fit[3] + 8*zeta(2)*zeta(3) + 3*zeta(5), 6))
