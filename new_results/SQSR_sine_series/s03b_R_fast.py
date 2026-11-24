# Faster high-precision T and R: Euler-Maclaurin with exact Taylor data (polygamma) at the cut T.
import sys, pickle, time
import sympy as sp
from mpmath import mp, mpf, psi, euler, zeta, log, pi, quad, inf, taylor, bernoulli, factorial, nstr
from symsum import t, LN2, GAM, ZETA, Hsym
DPS = int(sys.argv[1]); TT = int(sys.argv[2]); KK = int(sys.argv[3])
mp.dps = DPS + 30
law = pickle.load(open('b2_law.pkl', 'rb'))
H11, H12, H21, H22 = Hsym(1, 1), Hsym(1, 2), Hsym(2, 1), Hsym(2, 2)
Xs = [H11, H12, H21, H22]
def beta_t(e):
    n = 4*t + e
    return 2/((n+1)**2*(n-1)) if e == 1 else -2/((n-1)**2*(n+1))
def Cpi(e):
    n = 4*t + e
    return (-1)**((e+1)//2) * n*(2*t+1) / (2*(n**2-1)**2)
class ratfun:
    # rational function of t with exact Taylor coefficients from partial fractions
    def __init__(self, c):
        c = sp.apart(sp.together(c), t)
        self.terms = []
        for term in sp.Add.make_args(c):
            num, den = sp.fraction(sp.factor(term))
            if not den.has(t):
                assert sp.simplify(term) == 0 or not term.has(t), term
                assert sp.simplify(term) == 0, ('polynomial part', term)
                continue
            fl = sp.factor_list(den); lead = fl[0]; facs = [(f, e) for f, e in fl[1] if f.has(t)]
            const = sp.Integer(1)
            for f, e in fl[1]:
                if not f.has(t): const *= f**e
            assert len(facs) == 1; f, s_ = facs[0]
            p = sp.Poly(f, t).coeffs()[0]; a = sp.Rational(sp.expand(f - p*t), p)
            coef = sp.simplify(num/(lead*const*p**s_))
            self.terms.append((mpf(str(sp.N(coef, mp.dps+20))), s_, mpf(a.p)/a.q))
    def __call__(self, x):
        return mp.fsum(c/(x+a)**s for c, s, a in self.terms)
    def taylor(self, x, N):
        out = []
        for m in range(N+1):
            out.append(mp.fsum(c*(-1)**m*mp.binomial(s+m-1, m)/(x+a)**(s+m) for c, s, a in self.terms))
        return out
def pieces(expr):
    """expr = r0 + sum r_j X_j ; returns list of (lambdified r_j, kind)"""
    expr = sp.expand(expr.subs({LN2: sp.log(2), ZETA[2]: sp.pi**2/6}))
    out = []
    rest = expr
    for kind, X in enumerate(Xs, start=1):
        c = expr.coeff(X); rest = sp.expand(rest - c*X); c = sp.together(c)
        out.append((ratfun(c), kind))
    rest = sp.together(sp.expand(rest))
    assert not rest.free_symbols - {t}
    out.append((ratfun(rest), 0))
    return out
def Xval(kind, x):
    if kind == 0: return mpf(1)
    if kind == 1: return psi(0, x+1)+euler
    if kind == 2: return psi(0, 2*x+1)+euler
    if kind == 3: return zeta(2)-psi(1, x+1)
    return zeta(2)-psi(1, 2*x+1)
def Xtaylor(kind, x, N):
    if kind == 0: return [mpf(1)] + [mpf(0)]*N
    if kind == 1: return [psi(0, x+1)+euler] + [psi(k, x+1)/factorial(k) for k in range(1, N+1)]
    if kind == 2: return [psi(0, 2*x+1)+euler] + [2**k*psi(k, 2*x+1)/factorial(k) for k in range(1, N+1)]
    if kind == 3: return [zeta(2)-psi(1, x+1)] + [-psi(k+1, x+1)/factorial(k) for k in range(1, N+1)]
    return [zeta(2)-psi(1, 2*x+1)] + [-2**k*psi(k+1, 2*x+1)/factorial(k) for k in range(1, N+1)]
def em(pcs, t0, T, K):
    a = lambda x: sum(r(x)*Xval(k, x) for r, k in pcs)
    s = mp.fsum(a(mpf(j)) for j in range(t0, T))
    s += quad(a, [T, 2*T, 8*T, 64*T, inf]) + a(mpf(T))/2
    N = 2*K
    tay = [mpf(0)]*(N+1)
    for r, k in pcs:
        rt = r.taylor(mpf(T), N); xt = Xtaylor(k, mpf(T), N)
        for m in range(N+1):
            tay[m] += mp.fsum(rt[q]*xt[m-q] for q in range(m+1))
    for k in range(1, K+1):
        s -= bernoulli(2*k)/(2*k) * tay[2*k-1]       # B_{2k}/(2k)! * a^{(2k-1)} = B_{2k}/(2k) * tay[2k-1]
    return s
t_0 = time.time()
res = {}
for key in ('T', 'R'):
    tot = 0
    for e, t0 in ((1, 1), (3, 0)):
        full = sp.expand(beta_t(e)*law[e])
        if key == 'R':
            full = sp.expand(beta_t(e)*(law[e] - beta_t(e)*LN2 - Cpi(e)*6*ZETA[2]))
        tot += em(pieces(full), t0, TT, KK)
    res[key] = tot
l2 = log(2)
b1p = (l2-1)/2; b1pp = mpf(3)/8 - l2 + l2**2/4 + 5*pi**2/96 - pi**4/768
N3 = (87*zeta(3) + 16*l2**3 - 96*l2**2 + 72*l2 + 10*pi**2*l2 - 20*pi**2)/16
b1ppp = N3/2 - 3*b1p*b1pp - 3*res['T']
sq = pi**4/768 + 5*pi**2/96 - mpf(5)/8
chk = res['T'] - (l2*sq + pi**2*(1-zeta(3))/64 + res['R'])
mp.dps = DPS
print('dps', DPS, 'T', TT, 'K', KK, 'time %.1f s' % (time.time()-t_0))
print('T =', nstr(res['T'], DPS)); print('R =', nstr(res['R'], DPS)); print("b1''' =", nstr(b1ppp, DPS))
print('consistency T - (ln2*sum b_n^2 + pi^2(1-zeta3)/64 + R) =', nstr(chk, 3))
open('R_%d_%d.txt' % (DPS, TT), 'w').write('T %s\nR %s\nb1ppp %s\n' % (nstr(res['T'], DPS), nstr(res['R'], DPS), nstr(b1ppp, DPS)))
