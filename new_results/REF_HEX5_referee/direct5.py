"""Referee's own direct computation of P^{(k)}(1), k <= 5, for the weighted p-family of a
centrally symmetric convex polygon.  No project code is used.

K_p = { x : sum_j w_j |k_j <v_j,x>|^p <= 1 },  w_j k_j = c_j chosen so that K_1 is the polygon.
Kink line j passes through the corner V_j (and -V_j); l_j(x) = c_j |cross(V_j, x)|, sum_j l_j = 1 on dK_1.
Then the j-th term equals w_j^{-d} l_j^{1+d}, d = p-1.

Polar form: P(p) = int_0^{2pi} F^{-1/p} sqrt(1 + (F_theta/(p F))^2) dtheta, F = sum_j w_j^{-d} l_j(u)^{1+d}.
Differentiation under the integral sign: at each node the integrand is a truncated power series in d
(orders 0..NORD-1), so the k-th coefficient times k! integrates to P^{(k)}(1).
Each sector between consecutive corner directions is split at its midpoint; each half is parametrised by
the angle phi from its own corner, so the vanishing term |V_j| sin(phi) is computed with no cancellation.
Quadrature: own tanh-sinh on [0, mid] in the form phi = mid * sigmoid(pi*sinh(t)), all orders at once,
step halving until two successive levels agree.
"""
import mpmath as mp

NORD = 6  # orders 0..5


def jmul(a, b):
    n = len(a)
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1)) for k in range(n)]


def jexp(a):
    n = len(a)
    b = [mp.exp(a[0])] + [mp.mpf(0)] * (n - 1)
    for k in range(1, n):
        b[k] = mp.fsum(i * a[i] * b[k - i] for i in range(1, k + 1)) / k
    return b


def jlog(a):
    n = len(a)
    b = [mp.log(a[0])] + [mp.mpf(0)] * (n - 1)
    for k in range(1, n):
        b[k] = (a[k] - mp.fsum(i * b[i] * a[k - i] for i in range(1, k)) / k) / a[0]
    return b


def setup(V, w):
    """V: corners V_0..V_{m-1} counter-clockwise (the rest are -V_j); w: weight of kink line through V_j."""
    m = len(V)
    V = [(mp.mpf(a) if not isinstance(a, mp.mpf) else a, mp.mpf(b) if not isinstance(b, mp.mpf) else b) for a, b in V]
    cr = lambda P, Q: P[0] * Q[1] - P[1] * Q[0]
    M = mp.matrix(m, m)
    for q in range(m):
        for j in range(m):
            M[q, j] = abs(cr(V[j], V[q]))
    c = mp.lu_solve(M, mp.matrix([1] * m))
    c = [c[j] for j in range(m)]
    al = [mp.atan2(v[1], v[0]) for v in V]
    for i in range(1, m):
        while al[i] <= al[i - 1]:
            al[i] += 2 * mp.pi
    assert al[m - 1] < al[0] + mp.pi
    nrm = [mp.sqrt(v[0] ** 2 + v[1] ** 2) for v in V]
    # sanity: sum of terms is 1 at each corner
    for q in range(m):
        s = mp.fsum(c[j] * abs(cr(V[j], V[q])) for j in range(m))
        assert abs(s - 1) < mp.mpf(10) ** (-mp.mp.dps + 10)
    return dict(m=m, V=V, w=[mp.mpf(x) for x in w], c=c, al=al, nrm=nrm)


def jet_at(H, base_line, base_angle, s, phi, nord=NORD):
    """Integrand jet at theta = base_angle + s*phi; line base_line vanishes at phi = 0."""
    m = H['m']
    one = [mp.mpf(1)] + [mp.mpf(0)] * (nord - 1)
    F = [mp.mpf(0)] * nord
    Ft = [mp.mpf(0)] * nord
    for j in range(m):
        cj = H['c'][j] * H['nrm'][j]
        if j == base_line:
            sn = mp.sin(phi)
            ds = s * mp.cos(phi)
        else:
            dl = base_angle + s * phi - H['al'][j]
            sv = mp.sin(dl)
            sn = abs(sv)
            ds = mp.cos(dl) * (1 if sv > 0 else -1)
        l = cj * sn
        lp = cj * ds
        mu = mp.log(l) - mp.log(H['w'][j])
        e = [mp.mpf(1)]
        for k in range(1, nord):
            e.append(e[-1] * mu / k)
        # term = l e^{d mu};  d/dtheta term = (1+d) e^{d mu} l'
        for k in range(nord):
            F[k] += l * e[k]
            Ft[k] += lp * (e[k] + (e[k - 1] if k >= 1 else 0))
    inv1pd = [mp.mpf((-1) ** k) for k in range(nord)]
    lF = jlog(F)
    invF = jexp([-x for x in lF])
    g = jmul(jmul(Ft, invF), inv1pd)
    inner = [one[k] + v for k, v in enumerate(jmul(g, g))]
    sq = jexp([x / 2 for x in jlog(inner)])
    rho = jexp([-x for x in jmul(lF, inv1pd)])
    return jmul(rho, sq)


def ts_half(f, mid, nord, tol, maxlev=12):
    """int_0^mid f(phi) dphi for a vector-valued f, tanh-sinh with phi = mid*sigmoid(pi sinh t)."""
    # range of t: phi down to 10^-(2 dps)
    T = mp.mpf(1)
    while True:
        z = mp.pi * mp.sinh(-T)
        if z < -(2 * mp.mp.dps + 10) * mp.log(10):
            break
        T += mp.mpf(1) / 4
    cache = {}

    def node(t):
        key = t
        if key in cache:
            return cache[key]
        z = mp.pi * mp.sinh(t)
        if z < 0:
            ez = mp.exp(z)
            sig = ez / (1 + ez)
        else:
            sig = 1 / (1 + mp.exp(-z))
        wgt = mid * sig * (1 - sig) * mp.pi * mp.cosh(t)
        phi = mid * sig
        if phi == 0 or wgt == 0:
            v = [mp.mpf(0)] * nord
        else:
            fv = f(phi)
            v = [wgt * x for x in fv]
        cache[key] = v
        return v

    h = mp.mpf(1) / 2
    prev = None
    for lev in range(maxlev):
        n = int(T / h) + 1
        tot = [mp.mpf(0)] * nord
        for i in range(-n, n + 1):
            v = node(i * h)
            for k in range(nord):
                tot[k] += v[k]
        tot = [h * x for x in tot]
        if prev is not None:
            err = max(abs(tot[k] - prev[k]) for k in range(nord))
            if err < tol:
                return tot, err, lev
        prev = tot
        h /= 2
    return tot, err, lev


def sectors(H):
    m = H['m']
    al = H['al']
    out = []
    for q in range(m):
        a = al[q]
        if q < m - 1:
            b = al[q + 1]
            bl = q + 1
        else:
            b = al[0] + mp.pi
            bl = 0
        out.append((q, a, b, bl))
    return out


def edge_derivs(H, q, tol=None, nord=NORD):
    """Contribution of the arc over edge q (V_q -> V_{q+1}) to P^{(k)}(1), k = 0..nord-1, and an error estimate."""
    if tol is None:
        tol = mp.mpf(10) ** (-mp.mp.dps + 8)
    (q, a, b, bl) = sectors(H)[q]
    mid = (b - a) / 2
    r1, e1, l1 = ts_half(lambda phi: jet_at(H, q, a, 1, phi, nord), mid, nord, tol)
    r2, e2, l2 = ts_half(lambda phi: jet_at(H, bl, b, -1, phi, nord), mid, nord, tol)
    res = [(r1[k] + r2[k]) * mp.factorial(k) for k in range(nord)]
    return res, max(e1, e2) * mp.factorial(nord - 1), (l1, l2)


def all_derivs(H, nord=NORD, tol=None):
    """P^{(k)}(1) = 2 * sum over m edges (central symmetry)."""
    tot = [mp.mpf(0)] * nord
    per = []
    for q in range(H['m']):
        res, err, lev = edge_derivs(H, q, tol, nord)
        per.append((res, err, lev))
        for k in range(nord):
            tot[k] += 2 * res[k]
    return tot, per


def scalar_P(H, p):
    """P(p) by plain quadrature of the scalar integrand (for finite differences)."""
    m = H['m']
    tot = mp.mpf(0)
    for (q, a, b, bl) in sectors(H):
        mid = (b - a) / 2
        for (base_line, base, s) in ((q, a, 1), (bl, b, -1)):
            def f(phi):
                F = mp.mpf(0)
                Ft = mp.mpf(0)
                for j in range(m):
                    cj = H['c'][j] * H['nrm'][j]
                    if j == base_line:
                        sn = mp.sin(phi); ds = s * mp.cos(phi)
                    else:
                        dl = base + s * phi - H['al'][j]
                        sv = mp.sin(dl); sn = abs(sv); ds = mp.cos(dl) * (1 if sv > 0 else -1)
                    l = cj * sn; lp = cj * ds
                    wj = H['w'][j]
                    t = wj ** (1 - p) * l ** p
                    F += t
                    Ft += p * t * lp / l
                g = Ft / (p * F)
                return F ** (-1 / p) * mp.sqrt(1 + g * g)
            r, e, l = ts_half(lambda phi: [f(phi)], mid, 1, mp.mpf(10) ** (-mp.mp.dps + 6))
            tot += r[0]
    return 2 * tot
