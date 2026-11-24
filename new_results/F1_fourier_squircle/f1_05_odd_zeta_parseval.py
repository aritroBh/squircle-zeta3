"""
F1 step 5: explicit Fourier-series identities yielding the odd zeta values from the squircle.

u(theta) = d/dp log rho_p |_{p=2} = d/dp h_p |_{p=2} = -(cos^2 log|cos| + sin^2 log|sin|)/2
        = (log2/2 - 1/4) - sum_{n>=1} cos(4n theta) / (4n(4n^2-1))           (f1_02, PROVED)
<f,g> := (2/pi) int_0^{pi/2} f g dtheta.   beta_{2j}(theta) := sum_{n>=1} cos(4n theta)/n^{2j}
       = (-1)^{j+1} (2pi)^{2j} B_{2j}(2 theta/pi) / (2 (2j)!)   (Bernoulli polynomial, 0 <= theta <= pi/2).
CLAIMS (PROVED in result.md):
   zeta(3)    = 16 <u,1> + 8 <u, beta_2>  =  (4/pi) int_0^{pi/2} u(theta) [ (4 theta - pi)^2 + 8 - pi^2/3 ] dtheta
   zeta(2j+1) = 8 <u, beta_{2j} - 4 beta_{2j-2}>      (j >= 2)
Also: (i) Hurwitz's identity P^2 - 4 pi A = 2 pi^2 sum (16k^2-1) a_k^2 checked at p = 1.2, 1.5;
      (ii) diamond-end consistency of the uniform boundary-layer formula with D'(tau=0) = -8 pi.
"""
import mpmath as mp

mp.mp.dps = 60
pi = mp.pi; lam = mp.log(2)

def u(t):
    c, s = mp.cos(t), mp.sin(t)
    return -(c * c * mp.log(abs(c)) + s * s * mp.log(abs(s))) / 2

def beta(j, t):
    x = 2 * t / pi
    return (-1) ** (j + 1) * (2 * pi) ** (2 * j) * mp.bernpoly(2 * j, x) / (2 * mp.factorial(2 * j))

def ip(f, g):
    return 2 / pi * mp.quad(lambda t: f(t) * g(t), [0, pi / 8, pi / 4, 3 * pi / 8, pi / 2])

# sanity: beta_2 = ((4t-pi)^2 - pi^2/3)/4
t0 = mp.mpf('0.3')
print("beta_2(t0) - ((4t0-pi)^2-pi^2/3)/4 =", mp.nstr(beta(1, t0) - ((4 * t0 - pi) ** 2 - pi ** 2 / 3) / 4, 3))

z3 = 16 * ip(u, lambda t: 1) + 8 * ip(u, lambda t: beta(1, t))
z3b = 4 / pi * mp.quad(lambda t: u(t) * ((4 * t - pi) ** 2 + 8 - pi ** 2 / 3), [0, pi / 4, pi / 2])
print("zeta(3) via <u, 16 + 8 beta_2> :", mp.nstr(z3, 50), " diff", mp.nstr(z3 - mp.zeta(3), 3))
print("zeta(3) via single integral   :", mp.nstr(z3b, 50), " diff", mp.nstr(z3b - mp.zeta(3), 3))
for j in range(2, 6):
    val = 8 * ip(u, lambda t: beta(j, t) - 4 * beta(j - 1, t))
    print(f"zeta({2*j+1}) via 8<u, beta_{2*j} - 4 beta_{2*j-2}> = {mp.nstr(val, 40)}  diff {mp.nstr(val - mp.zeta(2 * j + 1), 3)}")

# the same with Aritro's own polynomial x^2 (x = 4 theta - pi in [-pi, pi]), as stated in result.md
val = 4 / pi * mp.quad(lambda t: (4 * t - pi) ** 2 * u(t), [0, pi / 4, pi / 2]) - pi ** 2 / 3 * (lam - mp.mpf(1) / 2) + 8 * lam - 4
print("zeta(3) = (4/pi) int (4t-pi)^2 u - (pi^2/3)(log2 - 1/2) + 8 log2 - 4 :", mp.nstr(val - mp.zeta(3), 3))
# and with x^4: <u, x^4> contains zeta(5) with coefficient -(2j)!/4 *(-1)^j ... check the general formula
def ip_poly(j):
    return 2 / pi * mp.quad(lambda t: (4 * t - pi) ** (2 * j) * u(t), [0, pi / 4, pi / 2])
def predicted(j):
    # <u, x^{2j}> = (lam/2 - 1/4) pi^{2j}/(2j+1) - (1/8) sum_i r_i Z_i,
    # r_i = (-1)^{i+1} 2 (2j)! pi^{2j-2i}/(2j-2i+1)!,  Z_i = 4^i (2 lam - 1) - sum_{l=1}^i 4^{i-l} zeta(2l+1)
    tot = (lam / 2 - mp.mpf(1) / 4) * pi ** (2 * j) / (2 * j + 1)
    for i in range(1, j + 1):
        r = (-1) ** (i + 1) * 2 * mp.factorial(2 * j) * pi ** (2 * j - 2 * i) / mp.factorial(2 * j - 2 * i + 1)
        Z = 4 ** i * (2 * lam - 1) - mp.fsum(4 ** (i - l) * mp.zeta(2 * l + 1) for l in range(1, i + 1))
        tot -= r * Z / 8
    return tot
for j in range(1, 5):
    print(f"<u, x^{2*j}>: quad {mp.nstr(ip_poly(j), 30)}  formula diff {mp.nstr(ip_poly(j) - predicted(j), 3)}")

# ---- (i) Hurwitz identity at p = 1.2 and 1.5 ----
mp.mp.dps = 30
def a_k(p, k):
    q = p / (p - 1)
    h = lambda t: (abs(mp.cos(t)) ** q + abs(mp.sin(t)) ** q) ** (1 / q)
    return 8 / pi * mp.quad(lambda t: h(t) * mp.cos(4 * k * t), mp.linspace(0, pi / 4, 9))
for p in [mp.mpf('1.5'), mp.mpf('1.2')]:
    P = pi * a_k(p, 0)
    A = 4 * mp.gamma(1 + 1 / p) ** 2 / mp.gamma(1 + 2 / p)
    K = 120
    ssum = mp.fsum((16 * k * k - 1) * a_k(p, k) ** 2 for k in range(1, K + 1))
    print(f"p={p}: P^2-4piA = {mp.nstr(P*P - 4*pi*A, 22)};  2pi^2 sum_(k<={K}) = {mp.nstr(2*pi**2*ssum, 22)};  diff {mp.nstr(P*P-4*pi*A-2*pi**2*ssum, 3)}")

# ---- (ii) diamond end: first-order change of the deficit from the boundary layer ----
# exact: D(tau) = P^2 - 4 pi A = (32 - 8 pi) - 8 pi tau + O(tau^2)  (P = sqrt2(4 + O(tau^2)), A = 2(1+tau) + O(tau^2))
# boundary layer: sum (16k^2-1)(a_k^2 - A_k^2) ~ -(8/pi) tau log2 + (4/pi) tau K,  K = int_0^inf (1/y - csch y)^2 dy
Kint = mp.quad(lambda y: (1 / y - 1 / mp.sinh(y)) ** 2, [0, 1, 10, mp.inf])
print(f"\nK = int_0^inf (1/y - csch y)^2 dy = {mp.nstr(Kint, 25)};  2 log2 - 1 = {mp.nstr(2*lam-1, 25)}; diff {mp.nstr(Kint-(2*lam-1),3)}")
print("=> boundary-layer prediction for the tau-coefficient of sum (16k^2-1)(a_k^2-A_k^2):",
      mp.nstr(-(8 / pi) * lam + (4 / pi) * Kint, 20), "   exact from P, A: -4/pi =", mp.nstr(-4 / pi, 20))
# eta-quotient asymptotics used: F(tau) = sum_k 1/(k sinh(2 pi k tau)) = pi/(12 tau) - log 2 + pi tau/6 + O(e^{-pi/(2 tau)})... check
for tau in [mp.mpf('0.05'), mp.mpf('0.02')]:
    F = mp.nsum(lambda k: 1 / (k * mp.sinh(2 * pi * k * tau)), [1, mp.inf])
    print(f"tau={tau}: F - (pi/(12 tau) - log2 + pi tau/6) = {mp.nstr(F - (pi/(12*tau) - lam + pi*tau/6), 3)}")
