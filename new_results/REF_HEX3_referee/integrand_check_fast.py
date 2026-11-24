"""Same derivation as integrand_check.py, but with the edge parameters set to random exact rationals
(L2 a rational square), so that sympy works with polynomials in t, X, Y, W only.  Each run is an exact
comparison with E3_fixed.pkl at that parameter point (Schwartz-Zippel style test of a polynomial identity)."""
import sympy as sp, pickle, random, time
t, X, Y, W, d = sp.symbols('t X Y W d')
P = pickle.load(open('<project>/research/08_theory/K1e_hexagon_closed/E3_fixed.pkl', 'rb'))
O = 4
def trunc(e):
    e = sp.expand(e); return sum(e.coeff(d, k)*d**k for k in range(O))
def mul(a_, b_): return trunc(a_*b_)
def ser_exp(z):
    out = 1; p = 1
    for k in range(1, O):
        p = mul(p, z); out += p/sp.factorial(k)
    return trunc(out)
def ser_log1p(z):
    out = 0; p = 1
    for k in range(1, O):
        p = mul(p, z); out += sp.Rational((-1)**(k+1), k)*p
    return trunc(out)
random.seed(11)
for trial in range(4):
    t0 = time.time()
    R = lambda a, b: sp.Rational(random.randint(a, b), random.randint(1, 9))
    B = sp.Rational(random.randint(1, 9), 10); r = R(1, 15)
    L = R(1, 30); L2 = L**2; q0 = R(-20, 20); V02 = (q0**2 + R(1, 40))/L2
    lA0, lA1, lB = R(-30, 5), R(-30, 5), R(-30, 5)
    A1 = 1 - B; A0 = 1 - B*r; kap = r - 1
    a = [A0*t, A1*(1-t), B*(1+kap*t)]; ap = [A0, -A1, B*kap]; u = [lA0 + X, lA1 + Y, lB + W]
    E = [ser_exp(d*uj) for uj in u]
    F = trunc(sum(aj*Ej for aj, Ej in zip(a, E)))
    S = trunc(sum(apj*Ej for apj, Ej in zip(ap, E)))
    assert sp.expand(F.coeff(d, 0)) == 1 and sp.expand(S.coeff(d, 0)) == 0
    logF = ser_log1p(F - 1); invF = ser_exp(-logF)
    Phi = mul(S, invF)
    inv1pd = sum((-d)**k for k in range(O))
    Fpow = ser_exp(mul(-2*logF, inv1pd))
    yyp = q0 + L2*t; yy = V02 + 2*q0*t + L2*t**2
    G = mul(Fpow, trunc(L2 - 2*yyp*Phi + yy*mul(Phi, Phi)))
    assert sp.expand(G.coeff(d, 0)) == L2
    z = trunc((G - L2)/L2)
    sq = ser_exp(ser_log1p(z)/2)
    mine = sp.expand(L*sq.coeff(d, 3))
    sy = {s: sp.Symbol(s) for s in ['A0','A1','B','kap','L2','q0','V02','lA0','lA1','lB']}
    sub = {sy['A0']: A0, sy['A1']: A1, sy['B']: B, sy['kap']: kap, sy['L2']: L2, sy['q0']: q0, sy['V02']: V02,
           sy['lA0']: lA0, sy['lA1']: lA1, sy['lB']: lB}
    theirs = sp.expand(sum(co.subs(sub, simultaneous=True)*t**n*X**aa*Y**bb*W**cc for (n, aa, bb, cc), co in P['nl'] + P['loc']))
    print(f'trial {trial}: B={B} r={r} L={L} q0={q0}: own - E3_fixed = {sp.expand(mine - theirs)}  ({time.time()-t0:.1f}s)', flush=True)
