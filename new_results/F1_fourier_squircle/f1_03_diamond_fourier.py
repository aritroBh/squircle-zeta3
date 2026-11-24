"""
F1 step 3 (main theorem): Fourier coefficients of the support function h_p near the diamond p = 1.

h_p(theta) = (|cos|^q + |sin|^q)^(1/q), tau = 1/q = 1 - 1/p.  a_k(p) = (4/pi) int_0^{pi/2} h_p cos(4k theta).
Exact representation (u = tan theta = e^{-tau s}):
   a_k(tau) - A_k = (8/pi) tau int_0^inf [(1+e^-s)^tau - 1] G_k(tau s) ds,
   A_k = a_k(diamond) = -(-1)^k 4 sqrt2/(pi(16k^2-1)),
   G_k(x) = cos(4k th) cos(th)^2 sin(th),  th = arctan(e^-x)  ( = (-1)^k cos(2k gd x) e^-x (1+e^-2x)^-3/2 ).
Asymptotic expansion (Watson-type, see result.md):  a_k - A_k ~ sum_{n>=2} c_n(k) tau^n with
   c_n(k) = (8/pi) sum_{j=1}^{n-1} (1/j!) I_{n-1-j, j} g_{k, n-1-j},
   I_{m,j} = int_0^inf s^m log^j(1+e^-s) ds  (weight m+j+1),   g_{k,m} = [x^m] G_k(x).
CLAIM: c_2(k) = (-1)^k sqrt2 zeta(2)/pi,  c_3(k) = (-1)^k sqrt2 zeta(3)/pi  for every k.
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
pi = mp.pi; s2 = mp.sqrt(2); z = mp.zeta

def A(k):
    return -(-1) ** k * 4 * s2 / (pi * (16 * k * k - 1))

def Gk(k, x):
    th = mp.atan(mp.e ** (-x))
    return mp.cos(4 * k * th) * mp.cos(th) ** 2 * mp.sin(th)

def a_minus_A(k, tau):
    f = lambda s: (mp.e ** (tau * mp.log1p(mp.e ** (-s))) - 1) * Gk(k, tau * s)
    return 8 / pi * tau * mp.quad(f, [0, 2, 8, 30, 120, mp.inf])

def a_direct(k, tau):
    q = 1 / tau
    h = lambda t: (abs(mp.cos(t)) ** q + abs(mp.sin(t)) ** q) ** tau
    pts = [0, pi / 8, pi / 4 - 6 * tau, pi / 4 - tau, pi / 4]
    return 8 / pi * mp.quad(lambda t: h(t) * mp.cos(4 * k * t), pts)

NMAX = 8
Icache = {}
def I(m, j):
    if (m, j) not in Icache:
        Icache[(m, j)] = mp.quad(lambda s: s ** m * mp.log1p(mp.e ** (-s)) ** j, [0, 2, 8, 30, 120, mp.inf])
    return Icache[(m, j)]

def gcoef(k):
    return mp.taylor(lambda x: Gk(k, x), 0, NMAX)

def cn(k, n, g):
    return 8 / pi * mp.fsum(I(n - 1 - j, j) * g[n - 1 - j] / mp.factorial(j) for j in range(1, n))

# sanity: the three weight-2/3 integrals
print("I_{0,1} - pi^2/12 =", mp.nstr(I(0, 1) - pi ** 2 / 12, 3))
print("I_{1,1} - 3 zeta(3)/4 =", mp.nstr(I(1, 1) - 3 * z(3) / 4, 3))
print("I_{0,2} - zeta(3)/4 =", mp.nstr(I(0, 2) - z(3) / 4, 3))

print("\n== direct quadrature vs exact s-representation (tau = 0.02) ==")
for k in [0, 1, 3]:
    t = mp.mpf('0.02')
    print(f"k={k}: a_direct - A - (s-int) = {mp.nstr(a_direct(k, t) - A(k) - a_minus_A(k, t), 3)}")

print("\n== Taylor coefficients c_n(k) (normalized: (-1)^k pi c_n / sqrt2) ==")
table = {}
for k in range(0, 6):
    g = gcoef(k)
    cs = {n: cn(k, n, g) for n in range(2, NMAX + 1)}
    table[k] = cs
    norm = lambda n: (-1) ** k * pi * cs[n] / s2
    print(f"k={k}: c2' - zeta2 = {mp.nstr(norm(2) - z(2), 3)};  c3' - zeta3 = {mp.nstr(norm(3) - z(3), 3)};"
          f"  c4' = {mp.nstr(norm(4), 20)};  c4'/zeta4 = {mp.nstr(norm(4) / z(4), 20)}")

print("\n== check the asymptotic expansion: R_N(tau) = (a_k - A_k - sum_{n<=N} c_n tau^n)/tau^(N+1) ==")
N = 6
for k in [0, 1, 2, 5]:
    line = []
    for t in [mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]:
        val = a_minus_A(k, t)
        R = (val - mp.fsum(table[k][n] * t ** n for n in range(2, N + 1))) / t ** (N + 1)
        line.append(mp.nstr(R, 10))
    print(f"k={k}: R_6 at tau=.02,.01,.005: {line};  c_7(k) = {mp.nstr(table[k][7], 10)}")

print("\n== tau^4 and tau^5 coefficients: PSLQ (pure weight test) ==")
for k in range(0, 4):
    x4 = (-1) ** k * pi * table[k][4] / s2
    r4 = mp.pslq([x4, z(4)], maxcoeff=10 ** 6, maxsteps=10 ** 6)
    x5 = (-1) ** k * pi * table[k][5] / s2
    r5 = mp.pslq([x5, z(5), z(2) * z(3)], maxcoeff=10 ** 6, maxsteps=10 ** 6)
    def resid(x, r, b):
        return abs(r[0] * x + mp.fsum(ri * bi for ri, bi in zip(r[1:], b))) if r else None
    print(f"k={k}: tau^4: {r4} resid {mp.nstr(resid(x4, r4, [z(4)]), 3)};   tau^5 [x, z5, z2z3]: {r5} resid {mp.nstr(resid(x5, r5, [z(5), z(2)*z(3)]), 3)}")
print("prediction tau^4: -(33/16) - 7k^2 (in units of zeta(4))")
for k in range(0, 4):
    x4 = (-1) ** k * pi * table[k][4] / s2
    print(f"   k={k}: x4/zeta4 + 33/16 + 7k^2 = {mp.nstr(x4 / z(4) + mp.mpf(33) / 16 + 7 * k * k, 3)}")

print("\n== individual weight-4/5 integrals (are they single-zeta on their own?) ==")
l2 = mp.log(2)
B4 = [z(4), mp.polylog(4, 0.5), l2 ** 4, z(2) * l2 ** 2, z(3) * l2]
for (m, j) in [(1, 2), (0, 3), (2, 1)]:
    r = mp.pslq([I(m, j)] + B4, maxcoeff=10 ** 6, maxsteps=10 ** 6)
    print(f"I_{m},{j} = {mp.nstr(I(m, j), 25)}  PSLQ [I, z4, Li4(1/2), l2^4, z2 l2^2, z3 l2]: {r}")

print("\n== perimeter (k=0): P = pi a_0 ; compare tau-series with D11's exact f_n ==")
Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, PI = sp.symbols('Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10 pi')
fn = {}
with open('research/08_theory/D11_rhombus_resurgence/step1_exact.txt') as fh:
    for line in fh:
        if line.startswith('f_'):
            lhs, rhs = line.split('=', 1)
            fn[int(lhs.strip()[2:])] = sp.sympify(rhs.strip(), locals={'Z2': Z2, 'Z3': Z3, 'Z4': Z4, 'Z5': Z5, 'Z6': Z6,
                                                                         'Z7': Z7, 'Z8': Z8, 'Z9': Z9, 'Z10': Z10, 'pi': PI})
wt = {Z2: 2, Z3: 3, Z4: 4, Z5: 5, Z6: 6, Z7: 7, Z8: 8, Z9: 9, Z10: 10, PI: 1}
def weights(expr):
    expr = sp.expand(expr)
    ws = set()
    for term in sp.Add.make_args(expr):
        w = 0
        for base, e in term.as_powers_dict().items():
            if base in wt:
                w += wt[base] * e
        ws.add(w)
    return ws
zval = {Z2: z(2), Z3: z(3), Z4: z(4), Z5: z(5), Z6: z(6), Z7: z(7), Z8: z(8), Z9: z(9), Z10: z(10), PI: pi}
print("tau-coefficient of P/(4 sqrt2) is e_n = (-1)^n (f_n + f_{n-1}); its set of weights, and pi*c_n(0) - 4sqrt2 e_n:")
for n in range(1, max(fn) + 1):
    e = sp.expand((-1) ** n * (fn[n] + fn[n - 1]))
    ws = weights(e)
    s_ = f"n={n:2d}: weights {sorted(ws)}  e_n = {sp.nsimplify(e) if n <= 5 else '...'}"
    if 2 <= n <= NMAX:
        ev = mp.mpf(sp.N(e.subs(zval), 70))
        s_ += f"   |pi c_n(0) - 4 sqrt2 e_n| = {mp.nstr(abs(pi * table[0][n] - 4 * s2 * ev), 3)}"
    print(s_)
# d-variable coefficients are mixed weight (for contrast)
print("for contrast, D11's f_n themselves have weights:", {n: sorted(weights(fn[n])) for n in range(2, 6)})

print("\n== uniform (boundary-layer) form: a_k ~ A_k + (-1)^k (sqrt2/pi)[1/(4k^2) - pi tau/(2k sinh(2 pi k tau))] ==")
for t in [mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]:
    worst = 0
    for k in [1, 2, 5, 10, 20, 50, 100, 200]:
        lead = (-1) ** k * s2 / pi * (mp.mpf(1) / (4 * k * k) - pi * t / (2 * k * mp.sinh(2 * pi * k * t)))
        err = a_minus_A(k, t) - lead
        worst = max(worst, abs(err))
    print(f"tau={mp.nstr(t,3)}: max_k |a_k - A_k - leading| = {mp.nstr(worst, 5)};  /tau^3 = {mp.nstr(worst / t ** 3, 5)}")
for k in [5, 50, 200]:
    t = mp.mpf('0.01')
    print(f"   tau=0.01, k={k}: a_k = {mp.nstr(A(k) + a_minus_A(k, t), 12)}, A_k = {mp.nstr(A(k), 12)}, "
          f"A_k+leading = {mp.nstr(A(k) + (-1)**k*s2/pi*(mp.mpf(1)/(4*k*k) - pi*t/(2*k*mp.sinh(2*pi*k*t))), 12)}")
