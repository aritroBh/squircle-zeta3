"""
VERIFIER V4: assorted independent checks.
 (a) zeta(3) single-integral formula and zeta(5) formula, tanh-sinh in the variable u = tan(theta), split at u=1 and mirrored
     (different quadrature and variable from the solver).
 (b) a_k(v) (v = d/dp log rho_p at p = 1): w-integral formula of result.md vs theta-quadrature vs the PSLQ closed forms;
     sign of the pi/4 term.
 (c) the classical integrals I_{0,3}, I_{1,2} quoted in the tau^4 proof.
 (d) Section 6 linear homotopy: P(t), A(t) for the Minkowski mix of diamond and square.
"""
import mpmath as mp
mp.mp.dps = 50
pi = mp.pi; lam = mp.log(2); G = mp.catalan; z3 = mp.zeta(3)

# (a) in u = tan th on [0,1] and the mirror (th -> pi/2 - th maps u -> 1/u; u(th) symmetric, x^2 term not).
def u_fun(th):
    c, s = mp.cos(th), mp.sin(th)
    return -(c * c * mp.log(abs(c)) + s * s * mp.log(abs(s))) / 2
def integrand_u(uu, weight):
    th = mp.atan(uu)
    return u_fun(th) * weight(th) / (1 + uu * uu)
w3 = lambda th: (4 * th - pi) ** 2 + 8 - pi ** 2 / 3
I = mp.quad(lambda x: integrand_u(x, w3), [0, 1], method='tanh-sinh') + \
    mp.quad(lambda x: integrand_u(1 / x, w3) / (x * x), [0, 1], method='tanh-sinh')
print("(a) (4/pi) int u [(4t-pi)^2 + 8 - pi^2/3] - zeta(3) =", mp.nstr(4 / pi * I - z3, 3))
def beta(j, t):
    return (-1) ** (j + 1) * (2 * pi) ** (2 * j) * mp.bernpoly(2 * j, 2 * t / pi) / (2 * mp.factorial(2 * j))
w5 = lambda th: beta(2, th) - 4 * beta(1, th)
I5 = mp.quad(lambda x: integrand_u(x, w5), [0, 1], method='tanh-sinh') + \
     mp.quad(lambda x: integrand_u(1 / x, w5) / (x * x), [0, 1], method='tanh-sinh')
print("    8<u, beta_4 - 4 beta_2> - zeta(5) =", mp.nstr(8 * 2 / pi * I5 - mp.zeta(5), 3))

# (b)
def v_theta(t):
    C, S = abs(mp.cos(t)), abs(mp.sin(t))
    return mp.log(C + S) - (C * mp.log(C) + S * mp.log(S)) / (C + S)
def ak_quad(k):
    return 4 / pi * mp.quad(lambda t: v_theta(t) * mp.cos(4 * k * t), [0, pi / 8, pi / 4, 3 * pi / 8, pi / 2])
def ak_w(k):
    f = lambda w: mp.log((1 + w) / (1 - w)) * mp.im((1 + 1j * w) ** (4 * k)) / (1 + w * w) ** (2 * k)
    return (-1) ** k / (2 * pi * k) * mp.quad(f, [-1, 0, 1])
closed = {1: pi / 4 - 1, 2: -pi / 4 + mp.mpf(4) / 3 - 2 / pi, 3: pi / 4 - mp.mpf(7) / 5 + 16 / (9 * pi),
          4: -pi / 4 + mp.mpf(152) / 105 - 98 / (45 * pi)}
for k in range(1, 5):
    q, w = ak_quad(k), ak_w(k)
    print(f"(b) k={k}: theta-quad {mp.nstr(q, 30)}  w-formula diff {mp.nstr(q - w, 3)}  closed-form diff {mp.nstr(q - closed[k], 3)}"
          f"   coefficient of pi/4 in closed form = {'+1' if k % 2 else '-1'}  ((-1)^k = {(-1)**k:+d})")
a0 = 4 / pi * mp.quad(v_theta, [0, pi / 4, pi / 2])
print("    k=0: a_0(v) - (4G/pi + log2 - pi/4) =", mp.nstr(a0 - (4 * G / pi + lam - pi / 4), 3))

# (c)
Li4h = mp.polylog(4, mp.mpf(1) / 2); z2, z4 = mp.zeta(2), mp.zeta(4)
I03 = mp.quad(lambda x: mp.log1p(x) ** 3 / x, [0, 1])
I12 = mp.quad(lambda s: s * mp.log1p(mp.e ** (-s)) ** 2, [0, 1, 10, 50, mp.inf])
I03c = 6 * z4 + mp.mpf(3) / 2 * z2 * lam ** 2 - mp.mpf(21) / 4 * z3 * lam - lam ** 4 / 4 - 6 * Li4h
I12c = -mp.mpf(15) / 4 * z4 + 4 * Li4h + lam ** 4 / 6 - z2 * lam ** 2 + mp.mpf(7) / 2 * z3 * lam
print("(c) I_{0,3} - closed =", mp.nstr(I03 - I03c, 3), ";  I_{1,2} - closed =", mp.nstr(I12 - I12c, 3))
# c4 bracket assembled from these exact constants (k = 0 part and k^2 part)
g0, g1 = mp.sqrt(2) / 4, mp.sqrt(2) / 8
for k in [0, 1, 4]:
    g2 = -(5 * mp.sqrt(2) / 32 + mp.sqrt(2) / 2 * k * k)
    c4 = 8 / pi * (g0 * I03c / 6 + g1 * I12c / 2 + g2 * mp.mpf(7) / 4 * z4)
    print(f"    k={k}: pi c4/sqrt2 + (33/16+7k^2) zeta4 = {mp.nstr(pi * c4 / mp.sqrt(2) + (mp.mpf(33)/16 + 7*k*k) * z4, 3)}")

# (d) Minkowski mix: h_t = (1-t) max(|c|,|s|) + t(|c|+|s|); P = int h, A = (1/2) int (h^2 - h'^2)
for t in [mp.mpf('0.3'), mp.mpf('0.7')]:
    h = lambda th: (1 - t) * max(abs(mp.cos(th)), abs(mp.sin(th))) + t * (abs(mp.cos(th)) + abs(mp.sin(th)))
    dh = lambda th: (1 - t) * (-mp.sin(th) if th < pi / 4 else mp.cos(th)) + t * (-mp.sin(th) + mp.cos(th))
    P = 8 * mp.quad(h, [0, pi / 4, pi / 2]) / 2
    A = 4 * mp.quad(lambda th: (h(th) ** 2 - dh(th) ** 2) / 2, [0, pi / 4, pi / 2])
    print(f"(d) t={t}: P - (4sqrt2(1-t)+8t) = {mp.nstr(P - (4*mp.sqrt(2)*(1-t) + 8*t), 3)};  A - (2+4t-2t^2) = {mp.nstr(A - (2 + 4*t - 2*t*t), 3)}")
