"""Shared definitions for the general-angle corner problem (Ahmed, Open Problem 2).

Notation (theta = turning angle, 0 < theta < pi):
  c = cos(theta), z = sin^2(theta/2) = (1-c)/2, w = cot^2(theta/2) = (1-z)/z, ell = log(1-z) = 2 log cos(theta/2)
  r(s) = sqrt(1 + 2 c s + s^2)
  h_c(s) = 1 + s - r(s)                       (Ahmed's defect density, eq. (hc))
  k_c(s) = r(s) - s - c = (1-c^2)/(s + c + r(s)) = (1-c) - h_c(s)   (cancellation-free form)
  M_theta(eps) = int_0^oo k_c(s) s^(eps-1) ds,   0 < eps < 1
  R_theta(eps) = eps * M_theta(eps)
"""
import mpmath as mp


def setdps(d):
    mp.mp.dps = d


def kc(s, c):
    r = mp.sqrt(1 + 2*c*s + s*s)
    return (1 - c*c)/(s + c + r)


def hc(s, c):
    r = mp.sqrt(1 + 2*c*s + s*s)
    return 2*(1 - c)*s/(1 + s + r)


# ---------------- closed forms of R_theta(eps) = eps*M_theta(eps) ----------------

def R_hyp(eps, th):
    """2 sin^2(th/2) * (pi eps / sin(pi eps)) * 2F1(-eps, 1+eps; 2; sin^2(th/2))"""
    z = mp.sin(th/2)**2
    return 2*z*(mp.pi*eps/mp.sin(mp.pi*eps))*mp.hyp2f1(-eps, 1 + eps, 2, z)


def R_ferrers(eps, th):
    """pi eps sin(th) P^{-1}_eps(cos th) / sin(pi eps)   (Ferrers function, mpmath legenp type=2)"""
    c = mp.cos(th)
    return mp.pi*eps*mp.sin(th)*mp.legenp(eps, -1, c, type=2)/mp.sin(mp.pi*eps)


def R_legdiff(eps, th):
    """pi eps (P_{eps-1}(c) - P_{eps+1}(c)) / ((2 eps + 1) sin(pi eps))"""
    c = mp.cos(th)
    return mp.pi*eps*(mp.legendre(eps - 1, c) - mp.legendre(eps + 1, c))/((2*eps + 1)*mp.sin(mp.pi*eps))


def R_deriv(eps, th):
    """pi sin^2(th) P'_eps(cos th) / ((1+eps) sin(pi eps)),  P' by numerical differentiation in x"""
    c = mp.cos(th)
    dP = mp.diff(lambda x: mp.legendre(eps, x), c)
    return mp.pi*mp.sin(th)**2*dP/((1 + eps)*mp.sin(mp.pi*eps))


def R_ahmed(eps):
    """Ahmed's R(eps) = Gamma(1+eps/2) Gamma((1-eps)/2) / (sqrt(pi) (1+eps))  (theta = pi/2)"""
    return mp.gamma(1 + eps/2)*mp.gamma((1 - eps)/2)/(mp.sqrt(mp.pi)*(1 + eps))


# ---------------- tail coefficients gamma_m ----------------

def gamma_m(m, c):
    """Coefficient of s^(1-m) in k_c(s) for s>1:  (-1)^m (P_{m-2}(c) - P_m(c)) / (2m-1),  m >= 2."""
    return (-1)**m*(mp.legendre(m - 2, c) - mp.legendre(m, c))/(2*m - 1)


def tail_series(n, th, tol=None):
    """S_theta(n) = sum_{m>=2} gamma_m eps n^(1-m)/(m-1-eps)
                  = sum_{m>=2} gamma_m / ( n^(m-1) ((m-1)(n-1) - 1) )."""
    c = mp.cos(th)
    eps = 1/(n - 1)
    if tol is None:
        tol = mp.mpf(10)**(-mp.mp.dps - 5)
    # Legendre polynomials by recurrence (exact, fast)
    P = [mp.mpf(1), c]
    S = mp.mpf(0)
    m = 2
    npow = n  # n^(m-1)
    while True:
        while len(P) <= m:
            k = len(P) - 1
            P.append(((2*k + 1)*c*P[k] - k*P[k - 1])/(k + 1))
        g = (-1)**m*(P[m - 2] - P[m])/(2*m - 1)
        term = g/(npow*((m - 1)*(n - 1) - 1))
        S += term
        # rigorous stop: |gamma_m| <= 2 for all m (Cauchy estimate), so the tail after m is
        # bounded by the geometric sum of 2 / (n^(m-1) |(m-1)(n-1)-1|)
        bound = 2/(npow*n*abs(m*(n - 1) - 1))*n/(n - 1)
        if m > 3 and m*(n - 1) > 2 and bound < tol:
            break
        m += 1
        npow *= n
        if m > 100000:
            raise RuntimeError("tail series did not converge")
    return S


def L_exact(n, th):
    """Exact formula: L_theta(n) = 1 + cos th + R_theta(eps) n^(-eps) - S_theta(n)."""
    eps = 1/(mp.mpf(n) - 1)
    return 1 + mp.cos(th) + R_hyp(eps, th)*mp.power(n, -eps) - tail_series(mp.mpf(n), th)


# ---------------- direct numerical lengths ----------------

def L_direct_x(n, th):
    """Literal arc length int_0^1 |u + n x^(n-1) v| dx, |u|=|v|=1, u.v = cos th.
    Breakpoints at x where the slope s = n x^(n-1) equals 2^j."""
    n = mp.mpf(n)
    c = mp.cos(th)
    eps = 1/(n - 1)

    def f(x):
        if x == 0:
            return mp.mpf(1)
        yp = n*mp.exp((n - 1)*mp.log(x))
        return mp.sqrt(1 + 2*c*yp + yp*yp)
    pts = [mp.mpf(0)]
    j = -60
    while True:
        s = mp.mpf(2)**j
        if s >= n:
            break
        x = (s/n)**eps
        if x > pts[-1]:
            pts.append(x)
        j += 1
    pts.append(mp.mpf(1))
    return mp.quad(f, pts)


def L_direct_s(n, th):
    """Defect identity (Ahmed Lemma 1, a change of variables x -> s = n x^(n-1)):
       2 - L = eps n^(-eps) int_0^n h_c(s) s^(eps-1) ds, evaluated with s = e^u."""
    n = mp.mpf(n)
    c = mp.cos(th)
    eps = 1/(n - 1)
    f = lambda u: hc(mp.exp(u), c)*mp.exp(eps*u)
    U = mp.log(n)
    pts = [-mp.inf, -40, -10, -3, 0]
    t = mp.mpf(3)
    while t < U:
        pts.append(t)
        t *= 2
    pts.append(U)
    I = mp.quad(f, pts)
    return 2 - eps*mp.power(n, -eps)*I


def digits(a, b):
    d = abs(a - b)
    if d == 0:
        return mp.inf
    return -mp.log10(d/abs(b)) if b != 0 else -mp.log10(d)
