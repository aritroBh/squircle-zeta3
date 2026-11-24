"""Taylor coefficients M_k(a) of the corner function
   M(u;a) = int_0^1 [ (1+g)^(-1-u) sqrt(1 + g^(2-2u) - 2cos(a) g^(1-u)) - 1 ] g^(u-1) dg = sum_k M_k(a) u^k
computed by expanding the integrand as a truncated power series in u at each quadrature node (exact
formal expansion, no numerical differentiation), then integrating each coefficient with tanh-sinh.
Own implementation; does not reuse the lens scripts.
"""
import mpmath as mp


def _mul(A, B, N):
    return [mp.fsum(A[i]*B[n - i] for i in range(n + 1)) for n in range(N + 1)]


def _sqrt(Q, N):
    b = [mp.sqrt(Q[0])]
    for n in range(1, N + 1):
        s = Q[n] - mp.fsum(b[k]*b[n - k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b


def integrand_series(g, c, N):
    L = mp.log(g)
    m = mp.log1p(g)
    fac = [mp.factorial(k) for k in range(N + 1)]
    E1 = [(-L)**k/fac[k] for k in range(N + 1)]        # g^(-u)
    E2 = [(-2*L)**k/fac[k] for k in range(N + 1)]      # g^(-2u)
    Q = [1 + g*g - 2*c*g] + [g*g*E2[k] - 2*c*g*E1[k] for k in range(1, N + 1)]
    B = _sqrt(Q, N)
    A = [(-m)**k/fac[k]/(1 + g) for k in range(N + 1)]  # (1+g)^(-1-u)
    AB = _mul(A, B, N)
    AB[0] -= 1
    G = [L**k/fac[k]/g for k in range(N + 1)]           # g^(u-1)
    return _mul(AB, G, N)


def corner_coeffs(a, N, pts=None):
    """Return [M_0(a), ..., M_N(a)] at the current mp precision."""
    c = mp.cos(a)
    cache = {}

    def F(g, k):
        key = g
        if key not in cache:
            cache[key] = integrand_series(g, c, N)
        return cache[key][k]
    if pts is None:
        pts = [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1]
    return [mp.quad(lambda g: F(g, k), pts) for k in range(N + 1)]


def M_exact(u, a, pts=None):
    c = mp.cos(a)
    if pts is None:
        pts = [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1]
    f = lambda g: ((1 + g)**(-1 - u)*mp.sqrt(1 + g**(2 - 2*u) - 2*c*g**(1 - u)) - 1)*g**(u - 1)
    return mp.quad(f, pts)
