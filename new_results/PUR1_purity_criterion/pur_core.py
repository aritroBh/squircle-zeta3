"""PUR1 core: symbolic and numeric machinery for the purity criterion.

F(tau) = int_0^V phi(v) [(1+e^{-2v/tau})^tau - 1] dv,  phi = e^v psi.
(a)  F_n = sum_k phi_k k! gamma_{k, n-1-k}
(b)  even part of psi:  tau * sum_k psi_{2k} tau^{2k} * (1/2) d_lambda^{2k} S(lambda,tau)|_{lambda=tau}
     odd part of psi :  tau * sum_k psi_{2k+1} tau^{2k+1} d_lambda^{2k+1} G(lambda,tau)|_{lambda=tau}
S(mu+x, mu) = nu (1-Q)/(nu^2-y^2), nu = mu/2, y = x/2,
log Q = sum_{k>=2} zeta(k)/k [(nu+y)^k + (nu-y)^k - (2nu)^k].
"""
import sympy as sp
import mpmath as mp

nu, y, t, v = sp.symbols('nu y tau v')
Z = {k: sp.Symbol('z%d' % k) for k in range(2, 30)}          # zeta(k) as symbols


def gam(a, b):
    return sp.Symbol('g_%d_%d' % (a, b))


def S_series(N):
    """S(nu,y) as polynomial, all homogeneous parts of degree <= N-1 (zeta-weight <= N)."""
    L = sum(Z[k] / k * ((nu + y)**k + (nu - y)**k - (2 * nu)**k) for k in range(2, N + 2))
    # exp(L) truncated at total degree N+1 in (nu,y)
    eps = sp.Symbol('eps')
    Le = sp.expand(L.subs({nu: eps * nu, y: eps * y}))
    Q = sp.series(sp.exp(Le), eps, 0, N + 2).removeO()
    Q = sp.expand(Q)
    res = 0
    for d in range(2, N + 2):
        comp = sp.expand(Q.coeff(eps, d))            # homogeneous degree d in (nu,y)
        q, r = sp.div(sp.Poly(-comp, nu, y), sp.Poly(nu**2 - y**2, nu, y))
        assert r.is_zero, "not divisible"
        res += nu * q.as_expr()
    return sp.expand(res)


def even_part_coeffs(psi_even, N, Sser=None):
    """psi_even: dict k->psi_{2k}. Returns dict n->coefficient of tau^n (n<=N) of the even-part sum."""
    if Sser is None:
        Sser = S_series(N)
    tot = 0
    for k, c in psi_even.items():
        if 2 * k + 1 > N:
            continue
        coef_y = Sser.coeff(y, 2 * k) if k > 0 else Sser.subs(y, 0)
        coef_y = sp.expand(coef_y)
        term = sp.Rational(1, 2) * sp.factorial(2 * k) * sp.Rational(1, 4**k) * coef_y.subs(nu, t / 2)
        tot += c * t**(2 * k + 1) * term
    tot = sp.expand(tot)
    return {n: sp.expand(tot.coeff(t, n)) for n in range(0, N + 1)}


def odd_part_coeffs(psi_odd, N):
    """psi_odd: dict k(odd)->psi_k. Coefficients (in gamma symbols) of tau * psi_k tau^k d^k G(tau,tau)."""
    out = {n: 0 for n in range(N + 1)}
    for k, c in psi_odd.items():
        # tau^{k+1} * sum_{a>=k,b>=1} gamma_ab a!/(a-k)! tau^{a-k+b}
        for a in range(k, N + 1):
            for b in range(1, N + 1):
                n = k + 1 + a - k + b
                if n <= N:
                    out[n] += c * sp.factorial(a) / sp.factorial(a - k) * gam(a, b)
    return {n: sp.expand(e) for n, e in out.items()}


def direct_coeffs(phi, N):
    """(a): F_n = sum_k phi_k k! gamma_{k,n-1-k}; phi dict k->phi_k."""
    out = {}
    for n in range(N + 1):
        e = 0
        for k in range(0, n - 1):
            b = n - 1 - k
            if b >= 1:
                e += phi.get(k, 0) * sp.factorial(k) * gam(k, b)
        out[n] = sp.expand(e)
    return out


def taylor(expr, var, n):
    s = sp.series(expr, var, 0, n + 1).removeO()
    return {k: sp.nsimplify(sp.expand(s).coeff(var, k)) for k in range(n + 1)}


# ---------- numerics ----------
_gcache = {}


def gamma_num(a, b):
    key = (a, b, mp.mp.dps)
    if key not in _gcache:
        f = lambda s: s**a * mp.log1p(mp.exp(-2 * s))**b
        pts = [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256]
        _gcache[key] = mp.quad(f, pts) / (mp.factorial(a) * mp.factorial(b))
    return _gcache[key]


def zeta_subs():
    return {Z[k]: mp.zeta(k) for k in range(2, 30)}


def num(expr):
    """evaluate an expression in z_k and g_a_b symbols (and pi, sqrt) numerically with mpmath."""
    syms = sorted(expr.free_symbols, key=lambda q: q.name)
    vals = []
    for s in syms:
        n = s.name
        if n.startswith('z'):
            vals.append(mp.zeta(int(n[1:])))
        elif n.startswith('g_'):
            _, a, b = n.split('_')
            vals.append(gamma_num(int(a), int(b)))
        else:
            raise ValueError(n)
    f = sp.lambdify(syms, expr, 'mpmath')
    return mp.mpf(f(*vals))


def to_pi(expr):
    """replace even zetas by rational multiples of pi^{2k}."""
    rep = {Z[2 * k]: sp.zeta(2 * k) for k in range(1, 14)}
    return sp.expand(expr.subs(rep))
