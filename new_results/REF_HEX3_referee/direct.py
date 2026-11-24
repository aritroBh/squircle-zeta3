"""Referee's own direct computation of P'''(1) for a weighted hexagon p-family.

K_p = { x : sum_j w_j |k_j <v_j, x>|^p <= 1 },  w_j k_j = c_j fixed by K_1 = given hexagon.
Equivalently F_p(x) = sum_j w_j^{1-p} |l_j(x)|^p with l_j(x) = c_j cross(V_j, x), sum_j |l_j| = 1 on dK_1.

Polar form:  P(p) = int_0^{2pi} F^{-1/p} sqrt(1 + (F_theta/(p F))^2) dtheta,  F = F_p(u_theta).
Method: "jet" = differentiate under the integral sign.  The integrand is evaluated as a truncated
power series in d = p-1 (order 3) at each node, and its d^3 coefficient (times 3!) is integrated with
tanh-sinh, split at every kink direction, each sector split in two halves integrated from its corner
so that the sine of the angle to the corner is computed without cancellation.
No code from the project is used.
"""
import mpmath as mp

N = 4  # jet length: orders 0..3

def jmul(a, b):
    return [mp.fsum(a[i] * b[n - i] for i in range(n + 1)) for n in range(N)]

def jexp(a):
    b = [mp.e ** a[0]] + [mp.mpf(0)] * (N - 1)
    for n in range(1, N):
        b[n] = mp.fsum(k * a[k] * b[n - k] for k in range(1, n + 1)) / n
    return b

def jlog(a):
    b = [mp.log(a[0])] + [mp.mpf(0)] * (N - 1)
    for n in range(1, N):
        b[n] = (a[n] - mp.fsum(k * b[k] * a[n - k] for k in range(1, n)) / n) / a[0]
    return b

def jscale(a, s):
    return [s * x for x in a]

def jadd(*args):
    return [mp.fsum(x[i] for x in args) for i in range(N)]

INV1PD = None  # 1/(1+d)

def setup(V, w):
    """V: three corners V0,V1,V2 (counter-clockwise, then -V0,...). w: weights of kink lines through V0,V1,V2."""
    V = [(mp.mpf(a), mp.mpf(b)) for a, b in V]
    cr = lambda P, Q: P[0] * Q[1] - P[1] * Q[0]
    M = mp.matrix(3, 3)
    for q in range(3):
        for j in range(3):
            M[q, j] = abs(cr(V[j], V[q]))
    c = mp.lu_solve(M, mp.matrix([1, 1, 1]))
    c = [c[0], c[1], c[2]]
    al = [mp.atan2(v[1], v[0]) for v in V]
    # make angles increasing
    for i in range(1, 3):
        while al[i] <= al[i - 1]:
            al[i] += 2 * mp.pi
    assert al[2] < al[0] + mp.pi
    norms = [mp.sqrt(v[0] ** 2 + v[1] ** 2) for v in V]
    return dict(V=V, w=[mp.mpf(x) for x in w], c=c, al=al, nrm=norms)

def integrand_jet(H, base_line, base_angle, s, phi):
    """theta = base_angle + s*phi; base_line is the index of the kink line through the endpoint corner."""
    one = [mp.mpf(1), mp.mpf(0), mp.mpf(0), mp.mpf(0)]
    terms, derivs = [], []
    for j in range(3):
        cj = H['c'][j] * H['nrm'][j]
        if j == base_line:
            # theta - al_j = s*phi + m*pi  -> |sin| = sin(phi), d|l|/dtheta = s*cos(phi)*sgn... (computed exactly)
            sn = mp.sin(phi)
            dabs = s * mp.cos(phi)  # d/dtheta of sin(phi) with phi = s*(theta-base)
        else:
            dlt = base_angle + s * phi - H['al'][j]
            sv = mp.sin(dlt)
            sn = abs(sv)
            dabs = mp.cos(dlt) * (1 if sv > 0 else -1)
        L = cj * sn
        Ld = cj * dabs
        m = mp.log(L) - mp.log(H['w'][j])
        e = [mp.mpf(1), m, m * m / 2, m ** 3 / 6]  # exp(d m)
        terms.append(jscale(e, L))
        derivs.append(jscale(jmul(e, [mp.mpf(1), mp.mpf(1), 0, 0]), Ld))  # (1+d) e^{dm} Ld
    F = jadd(*terms)
    Ft = jadd(*derivs)
    inv1pd = [mp.mpf(1), mp.mpf(-1), mp.mpf(1), mp.mpf(-1)]
    lF = jlog(F)
    invF = jexp(jscale(lF, -1))
    g = jmul(jmul(Ft, invF), inv1pd)
    inner = jadd(one, jmul(g, g))
    sq = jexp(jscale(jlog(inner), mp.mpf(1) / 2))
    rho = jexp(jscale(jmul(lF, inv1pd), -1))
    return jmul(rho, sq)

def edge_jet(H, q, order=3, maxdeg=None):
    """Integral over the sector of edge q (V_q -> V_{q+1}) of the d^order coefficient; returns list of all orders."""
    al = H['al']
    a = al[q]
    if q < 2:
        b = al[q + 1]; bl = q + 1
    else:
        b = al[0] + mp.pi; bl = 0
    mid = (b - a) / 2
    out = []
    for k in range(N):
        f1 = lambda phi, k=k: integrand_jet(H, q, a, 1, phi)[k]
        f2 = lambda phi, k=k: integrand_jet(H, bl, b, -1, phi)[k]
        kw = {} if maxdeg is None else {'maxdegree': maxdeg}
        I = mp.quad(f1, [0, mid], **kw) + mp.quad(f2, [0, mid], **kw)
        out.append(I * mp.factorial(k))  # derivative contribution
    return out

def edge_jet_fast(H, q, maxdeg=None):
    """Same, all 4 orders from a single quadrature pass (vectorised by caching nodes)."""
    al = H['al']
    a = al[q]
    if q < 2:
        b = al[q + 1]; bl = q + 1
    else:
        b = al[0] + mp.pi; bl = 0
    mid = (b - a) / 2
    res = [mp.mpf(0)] * N
    for (base_line, base, s) in ((q, a, 1), (bl, b, -1)):
        cache = {}
        def f(phi, k):
            key = phi
            if key not in cache:
                cache[key] = integrand_jet(H, base_line, base, s, phi)
            return cache[key][k]
        kw = {} if maxdeg is None else {'maxdegree': maxdeg}
        for k in range(N):
            res[k] += mp.quad(lambda phi: f(phi, k), [0, mid], **kw)
    return [res[k] * mp.factorial(k) for k in range(N)]

def scalar_integrand(H, base_line, base_angle, s, phi, p):
    F = 0; Ft = 0
    for j in range(3):
        cj = H['c'][j] * H['nrm'][j]
        if j == base_line:
            sn = mp.sin(phi); dabs = s * mp.cos(phi)
        else:
            dlt = base_angle + s * phi - H['al'][j]
            sv = mp.sin(dlt); sn = abs(sv); dabs = mp.cos(dlt) * (1 if sv > 0 else -1)
        L = cj * sn; Ld = cj * dabs
        w = H['w'][j]
        F += w ** (1 - p) * L ** p
        Ft += w ** (1 - p) * p * L ** (p - 1) * Ld
    g = Ft / (p * F)
    return F ** (-1 / p) * mp.sqrt(1 + g * g)

def perimeter(H, p):
    al = H['al']; tot = 0
    for q in range(3):
        a = al[q]
        if q < 2:
            b = al[q + 1]; bl = q + 1
        else:
            b = al[0] + mp.pi; bl = 0
        mid = (b - a) / 2
        tot += mp.quad(lambda phi: scalar_integrand(H, q, a, 1, phi, p), [0, mid])
        tot += mp.quad(lambda phi: scalar_integrand(H, bl, b, -1, phi, p), [0, mid])
    return 2 * tot
