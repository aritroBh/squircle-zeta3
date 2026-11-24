# rf4: numerical checks of Theorem A(i), Lemma 3 (branch identity), Lemma 4 (Hankel bound, complex h),
# and the Gamma-ratio bound of Section 2.4 (own code).
import mpmath as mp
from flint import acb, arb, ctx
from rf_series import mul, log, exp
mp.mp.dps = 40

# ---- (a) Theorem A(i): Ehat = h sum_n C(h,n) I_n,  I_n = int_0^1 x^{a-1}(1+x)^{-1-h} dx,  a = n(1-h)/h + h ----
print("(a) Theorem A(i), Laplace form vs s-integral")
for p in (3, 4, 10):
    h = mp.mpf(1)/p; kap = h/(1-h)
    f = lambda s: mp.exp(-h*kap*s)*(1+mp.exp(-kap*s))**(-1-h)*mp.expm1(h*mp.log1p(mp.exp(-s)))
    direct = h*kap*mp.quad(f, [0, 1, 5, 20, 60, 150, 400, mp.inf])
    def term(n):
        n = int(n); a = n*(1-h)/h + h
        return mp.binomial(h, n)*mp.hyp2f1(1+h, a, a+1, -1)/a
    lap = h*mp.nsum(term, [1, mp.inf])
    print("  p=%d  direct=%s  laplace=%s  diff=%s" % (p, mp.nstr(direct, 20), mp.nstr(lap, 20), mp.nstr(direct-lap, 3)))

# ---- (b) Lemma 3: A(y,h) = (y0-y)^{-1-h} Phi(y-y0,h) on V_+, A continued analytically from y=0 ----
print("(b) Lemma 3 at complex h; ell(y) by integrating e^w/(1+e^w) along [0,y]")
y0 = mp.mpc(0, mp.pi); Ly0 = mp.log(y0)
def ell(y):
    g = lambda w: mp.exp(w)/(1+mp.exp(w))
    return mp.log(2) + mp.quad(lambda t: y*g(t*y), mp.linspace(0, 1, 41))
def A(y, h): return mp.exp(2*y - (1+h)*ell(y))
def rhs(y, h):
    u = y - y0
    pw = mp.exp((-1-h)*Ly0)*mp.exp((-1-h)*mp.log(1-y/y0))
    return pw*mp.exp(2*u - (1+h)*mp.log((mp.exp(u)-1)/u))
mx = 0
for h in (mp.mpc(0.2, 0.15), mp.mpc(-0.25, 0), mp.mpc(0.1, -0.3), mp.mpf(0.7)):
    for y in (mp.mpc(0.3, mp.pi), mp.mpc(-0.3, mp.pi), mp.mpc(0.01, 1.2*mp.pi), mp.mpc(-0.01, 1.2*mp.pi),
              mp.mpc(0.4, 0.8*mp.pi), mp.mpc(-1.0, 1.1*mp.pi), mp.mpc(1.2, 1.05*mp.pi), mp.mpc(0, 0.6*mp.pi)):
        if abs(y-y0) >= mp.pi/2: continue
        e = abs(A(y, h)/rhs(y, h) - 1); mx = max(mx, e); print("   h=%s y=%s  rel.err %s" % (mp.nstr(h, 3), mp.nstr(y, 4), mp.nstr(e, 3)))
print("  max |A/RHS - 1| over test points:", mp.nstr(mx, 3))

# ---- (c) Lemma 4: |[y^m] R_L(.,h)| pi^m m^{L-|h|} bounded, h on |h| = 1/4 ----
print("(c) Lemma 4: sup_m |[y^m]R_L| pi^m m^{L-|h|} over windows of m (should stay bounded)")
ctx.prec = 1200
NM = 401
def Acoef(hc):
    fac = [arb(1)]
    for n in range(1, NM+2): fac.append(fac[-1]*n)
    ey = [acb(1)/fac[n] for n in range(NM)]
    u = [ey[0]+1] + ey[1:]
    l = log(u, NM, acb(arb(2).log()))
    ex = [acb(2)*0 for _ in range(NM)]
    ex = [-(1+hc)*c for c in l]
    ex[1] += 2
    return exp(ex, NM)
def phil(l, h):   # phi_l(h) by Cauchy integral of Phi(u,h)
    f = lambda t: mp.exp(2*t-(1+h)*mp.log((mp.exp(t)-1)/t)) if t != 0 else 1
    return mp.quad(lambda th: f(mp.mpf(1.5)*mp.expj(th))*mp.expj(-l*th), [0, 2*mp.pi])/(2*mp.pi*mp.mpf(1.5)**l)
for theta in (0, mp.pi/2, mp.pi, 3*mp.pi/2, mp.pi/4):
    hm = 0.25*mp.expj(theta)
    hc = acb(mp.nstr(hm.real, 30), mp.nstr(hm.imag, 30))
    Ac = Acoef(hc)
    A_m = [mp.mpc(mp.mpf(c.real.mid().str(60, radius=False)), mp.mpf(c.imag.mid().str(60, radius=False))) for c in Ac]

    phv = [phil(l, hm) for l in range(3)]
    out = []
    for L in (1, 2, 3):
        mxw = {}
        for m in range(1, NM):
            Sm = 0
            for l in range(L):
                for yy0, LL in ((y0, Ly0), (-y0, mp.conj(Ly0))):
                    Sm += (-1)**l*phv[l]*yy0**(l-1-m)*mp.exp(-hm*LL)*mp.rf(1+hm-l, m)/mp.factorial(m)
            R = A_m[m] - Sm
            sc = abs(R)*mp.pi**m*mp.mpf(m)**(L-abs(hm))
            w = (m-1)//100
            mxw[w] = max(mxw.get(w, 0), sc)
        out.append("L=%d: " % L + " ".join(mp.nstr(mxw[w], 3) for w in sorted(mxw)))
    print("  h=0.25e^{i%.2f} | " % float(theta) + " | ".join(out))

# ---- (d) Gamma bound |Gamma(M+h)/Gamma(M)| <= 1.4 M^{|h|}, |h|<=1/4 ----
mxr = 0
for M in list(range(1, 60)) + [100, 1000, 10**5]:
    for j in range(64):
        hh = 0.25*mp.expj(2*mp.pi*j/64)
        mxr = max(mxr, abs(mp.gamma(M+hh)/mp.gamma(M))/mp.mpf(M)**0.25)
print("(d) max |Gamma(M+h)/Gamma(M)|/M^{1/4} on |h|=1/4 grid:", mp.nstr(mxr, 5), "(claim <= 1.4)")
