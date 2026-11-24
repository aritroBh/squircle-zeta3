import sympy as sp, tow2lib as T
d = sp.Symbol('d', real=True); m = sp.Symbol('m', real=True); s = sp.Symbol('s', positive=True, integer=True)
i = sp.I
# (1) telescoping with symbolic s via Gamma:  Pi_s(d)=2^s Gamma(d-1/2)/Gamma(d-s-1/2)
Pi = lambda dd: 2**s*sp.gamma(dd-sp.Rational(1,2))/sp.gamma(dd-s-sp.Rational(1,2))
ratio = sp.gammasimp(Pi(d+1)/Pi(d))
print('Pi_s(d+1)/Pi_s(d) (symbolic s) =', sp.simplify(ratio), ' ; minus (2d-1)/(2d-2s-1):', sp.simplify(ratio - (2*d-1)/(2*d-2*s-1)))
# (2) explicit e(2s,s), e(2s+1,s) top equations for s=1..10, symbolic d,m, using T.* definitions (tau from DEFINITION, not table)
def Pis(ss, dd):
    out = 1
    for j in range(1, ss+1): out *= (2*dd-2*j-1)
    return out
def etop(ss, dd, odd, Cn):
    c = (-1)**ss*Cn/Pis(ss, dd)
    return c if not odd else -i*(dd-ss-2)*c
A = T.A
okall = True
for ss in range(1, 11):
    C = A[ss+1]
    # even: Lambda^(2s)(s) = sum over parents delta'=2s only
    cA = T.tau1(0)*(m+1)/(2*d-2*ss-1); cB = T.tau2(0)*(d-m+1)/(2*d-2*ss-1)
    lam = cA*etop(ss, d+1, False, C)*T.rA0(d, m) + cB*etop(ss, d+1, False, C)*T.rB0(d, m)
    r1 = sp.simplify(sp.radsimp(sp.together(lam - etop(ss, d, False, C))))
    # odd: delta'=2s+1 (tau_0) and delta'=2s (tau_1)
    lam2 = 0
    for dp, odd in ((2*ss+1, True), (2*ss, False)):
        cA = T.tau1(2*ss+1-dp)*(m+1)/(2*d-dp-1); cB = T.tau2(2*ss+1-dp)*(d-m+1)/(2*d-dp-1)
        lam2 += cA*etop(ss, d+1, odd, C)*T.rA0(d, m) + cB*etop(ss, d+1, odd, C)*T.rB0(d, m)
    r2 = sp.simplify(sp.radsimp(sp.together(lam2 - etop(ss, d, True, C))))
    print(f's={ss}: even top remainder {r1}, odd top remainder {r2}')
    okall &= (r1 == 0 and r2 == 0)
print('all remainders zero for s=1..10 (symbolic d,m; tau from definitions):', okall)
# (3) conj via subs(I->-I) equals sp.conjugate for real d,m
for dl in (1, 2, 3, 4):
    a = T.arho_simple(dl)
    a_sub = a.subs(sp.I, -sp.I)
    # re-symbolise with real symbols
    a_r = a.subs({T.d: d, T.m: m})
    print('delta', dl, 'subs(I->-I) == conjugate:', sp.simplify(sp.expand(a_r.subs(sp.I, -sp.I) - sp.conjugate(a_r))) == 0)
for n in range(0, 7):
    print('conj(tau1_%d) == (-1)^n tau2_%d :' % (n, n), sp.simplify(sp.conjugate(T.tau1f(n)) - (-1)**n*T.tau2f(n)) == 0)
