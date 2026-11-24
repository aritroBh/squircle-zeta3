"""(1) Symbolic check of A2 = the '-w' part of CX6 (10) (sympy, exact Laplace integrals).
(2) Numerical check at sample s on Gamma: A_k(s) from CX6 (10) (full kernel), A3_k(s) from the referee's
    remainder kernel, and A_k - A3_k - A2_k(closed form) = 0; plus check of (B3) at those s.
(3) Direct U_k = -J_k + (2/pi) Im{I1 + I2 + I3 - int_0^{1/2} H E Q_k - int_0^{1/2} H E A2_k} at k = 100,
    computed from the definitions (I2 via the remainder kernel, no closed form for A3), to be compared
    with rest - S_k from r2_rest_table (which uses data M_k, quadrature M_main and the closed form A2)."""
import sympy as sp, mpmath as mp, sys
# ---- (1) symbolic
s, L = sp.symbols('s L', positive=True); k = sp.Symbol('k', positive=True, integer=True)
kk = 7
expr = (-(s*L)**kk - kk*(s*L)**(kk-1))*sp.exp(-(s+2)*L)/sp.factorial(kk)
val = sp.simplify(sp.integrate(expr, (L, 0, sp.oo)))
A2 = -2*(s+1)*s**(kk-1)/(s+2)**(kk+1)
print("(1) k=7 symbolic: integral of '-w' part minus closed form A2 =", sp.simplify(val - A2))
# ---- (2) numeric kernels
mp.mp.dps = 60
om = mp.expjpi(mp.mpf(1)/4)
def kernel_full(s, Lc, k):
    w = mp.exp(-Lc); b = mp.log(1+w)
    return mp.exp(-(s+1)*Lc)*((s*Lc-b)**k/(1+w) - (s*Lc)**k)
def kernel_rem(s, Lc, k):   # full kernel minus the '-w' parts that make A2
    w = mp.exp(-Lc); b = mp.log(1+w)
    return mp.exp(-(s+1)*Lc)*((s*Lc-b)**k/(1+w) - (s*Lc)**k + (s*Lc)**k*w + k*(s*Lc)**(k-1)*w)
def Lint(kern, s, k, shift):
    ph = mp.arg(s+shift); d = mp.expj(-ph)            # rotate L to arg = -arg(s+shift), inside |arg L| <= pi/4
    Lp = mp.mpf(k)/abs(s+shift)
    pts = [0, Lp/4, Lp/2, 3*Lp/4, Lp, 1.25*Lp, 1.5*Lp, 2*Lp, 3*Lp, 5*Lp, 8*Lp, 12*Lp]   # tail beyond 12Lp < e^{-5k}
    return mp.quad(lambda t: kern(s, t*d, k)*d, pts, maxdegree=6)/mp.factorial(k)
def A2c(s, k): return -2*(s+1)*s**(k-1)/(s+2)**(k+1)
def B3(s, k):
    rho = abs(s)
    v = 2*rho**k/abs(s+3)**(k+1) + 4*rho**(k-1)/abs(s+3)**k
    v += 2*mp.fsum(mp.mpf(2)**i/mp.factorial(i)*rho**(k-i)/abs(s+i+1)**(k-i+1) for i in range(2, k+1))
    return v
for k in (20, 40):
    for r in (0, 2, 10):
        sv = mp.mpf(1)/2 + r*om
        Af = Lint(kernel_full, sv, k, 2); A3 = Lint(kernel_rem, sv, k, 3)
        print(f"(2) k={k} r={r}: |A_k|={mp.nstr(abs(Af),6)} |A3|={mp.nstr(abs(A3),6)} "
              f"|A_k-A3-A2|/|A2|={mp.nstr(abs(Af-A3-A2c(sv,k))/abs(A2c(sv,k)),3)}  |A3|/B3={mp.nstr(abs(A3)/B3(sv,k),4)}", flush=True)
# ---- (3) direct U_k at k=100
k = int(sys.argv[1]) if len(sys.argv) > 1 else 100
N = k+1
mp.mp.dps = 40   # absolute accuracy suffices: integrand magnitudes <= 1, target 1e-25
pi = mp.pi
H = lambda s: mp.beta(s/2, mp.mpf(3)/2)/2
E = lambda s: mp.expj(pi*s/2)
K = lambda s: 1/(2*mp.cos(pi*s/2))
Q = lambda s: s**k/(1+s)**(k+1)
G = lambda r: mp.mpf(1)/2 + r*om
rpts = [0, 1, 2, 4, 8, 12, 16, 20, 30, 45, 70, 110]
I1 = mp.quad(lambda r: H(G(r))*(K(G(r))-E(G(r)))*A2c(G(r), k)*om, rpts, maxdegree=6)
I3 = mp.quad(lambda r: H(G(r))*(K(G(r))-E(G(r)))*Q(G(r))*om, rpts, maxdegree=6)
segQ = mp.quad(lambda x: H(x)*E(x)*Q(x), [0, mp.mpf(1)/4, mp.mpf(1)/2])
segA = mp.quad(lambda x: H(x)*E(x)*A2c(x, k), [0, mp.mpf(1)/4, mp.mpf(1)/2])
J = (-1)**k*mp.log(2)**N/mp.factorial(N)
print(f"(3) k={k}: I1={mp.nstr(I1,12)}  I3={mp.nstr(I3,12)}  segQ={mp.nstr(segQ,8)} segA={mp.nstr(segA,8)}", flush=True)
cnt = [0]
def f2(r):
    sv = G(r); cnt[0] += 1
    return H(sv)*K(sv)*Lint(kernel_rem, sv, k, 3)*om
I2 = mp.quad(f2, rpts, maxdegree=6)
print(f"    I2={mp.nstr(I2,15)}  (outer evals {cnt[0]})")
U = -J + 2/pi*mp.im(I1 + I2 + I3 - segQ - segA)
print(f"    U_k(direct) = {mp.nstr(U, 20)}")
print(f"    pieces (2/pi)Im: I1={mp.nstr(2/pi*mp.im(I1),8)} I2={mp.nstr(2/pi*mp.im(I2),8)} I3={mp.nstr(2/pi*mp.im(I3),8)}")
