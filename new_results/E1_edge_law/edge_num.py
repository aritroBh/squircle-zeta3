# edge_num.py: numerical tools for K_p = {sum_j w_j e^{d lam_j} |<n_j,x>|^p <= 1}, d = p-1.
# Kink j: n_j = (-sin b_j, cos b_j), so <n_j,(cos th, sin th)> = sin(th - b_j); kink line = direction b_j.
# Problem convention: lam_j = 0 (weights fixed).  "Power" convention sum |<v_j,x>|^p: w_j=|v_j|, lam_j = log|v_j|.
import mpmath as mp

def polygon(bs, ws):
    """vertices of K_1 in angular order on [0, 2pi): list of (theta, V, kink index)"""
    m = len(bs)
    def S0(th):
        return mp.fsum([w*abs(mp.sin(th - b)) for b, w in zip(bs, ws)])
    dirs = []
    for j, b in enumerate(bs):
        for k in (0, 1):
            th = (b + k*mp.pi) % (2*mp.pi)
            dirs.append((th, j))
    dirs.sort(key=lambda x: x[0])
    out = []
    for th, j in dirs:
        r = 1/S0(th)
        out.append((th, mp.matrix([r*mp.cos(th), r*mp.sin(th)]), j))
    return out

def edges(bs, ws, lams):
    """edge data for every edge of K_1 (counter-clockwise), as dicts"""
    vs = polygon(bs, ws)
    n = len(vs)
    E = []
    for k in range(n):
        th0, V0, j0 = vs[k]
        th1, V1, j1 = vs[(k+1) % n]
        D = V1 - V0
        L = mp.sqrt(D[0]**2 + D[1]**2)
        tau = D/L
        ne = mp.matrix([tau[1], -tau[0]])           # outward normal (ccw polygon)
        h = V0[0]*ne[0] + V0[1]*ne[1]
        s0 = V0[0]*tau[0] + V0[1]*tau[1]
        s1 = V1[0]*tau[0] + V1[1]*tau[1]
        nrm = [mp.matrix([-mp.sin(b), mp.cos(b)]) for b in bs]
        c0 = [nn[0]*V0[0] + nn[1]*V0[1] for nn in nrm]
        c1 = [nn[0]*V1[0] + nn[1]*V1[1] for nn in nrm]
        E.append(dict(V0=V0, V1=V1, j0=j0, j1=j1, L=L, h=h, s0=s0, s1=s1, c0=c0, c1=c1,
                      ws=ws, lams=lams, th0=th0, th1=th1,
                      alpha=ws[j0]*abs(c1[j0]), beta=ws[j1]*abs(c0[j1])))
    return E

# truncated power series in d (lists of length N)
N = 4
def smul(a, b): return [mp.fsum([a[i]*b[k-i] for i in range(k+1)]) for k in range(N)]
def sinv(a):
    r = [1/a[0]] + [mp.mpf(0)]*(N-1)
    for k in range(1, N):
        r[k] = -mp.fsum([a[i]*r[k-i] for i in range(1, k+1)])/a[0]
    return r
def slog1(a):
    x = [mp.mpf(0)] + a[1:]; res = [mp.mpf(0)]*N; pw = [mp.mpf(1)] + [mp.mpf(0)]*(N-1)
    for j in range(1, N):
        pw = smul(pw, x); res = [res[k] + mp.mpf((-1)**(j+1))/j*pw[k] for k in range(N)]
    return res
def sexp0(a):
    res = [mp.mpf(1)] + [mp.mpf(0)]*(N-1); pw = list(res)
    for j in range(1, N):
        pw = smul(pw, a); res = [res[k] + pw[k]/mp.factorial(j) for k in range(N)]
    return res

def E_series(e, t):
    """[E0, E1, E2, E3](t) for edge e: coefficients of d^k in |x_p'(t)|"""
    ws, lams, c0, c1 = e['ws'], e['lams'], e['c0'], e['c1']
    a, ap, l = [], [], []
    for i in range(len(ws)):
        v = (1-t)*c0[i] + t*c1[i]
        sg = mp.sign(v) if v != 0 else mp.sign(c0[i] + c1[i])
        a.append(ws[i]*abs(v)); ap.append(ws[i]*sg*(c1[i] - c0[i]))
        l.append(mp.log(abs(v)) + lams[i] if v != 0 else mp.mpf(0))
    F = [mp.fsum([a[i]*l[i]**k for i in range(len(ws))])/mp.factorial(k) for k in range(N)]
    Num = [mp.fsum([ap[i]*l[i]**k for i in range(len(ws))])/mp.factorial(k) for k in range(N)]
    Num[0] = mp.mpf(0)
    Phi = smul(Num, sinv(F))
    lF = slog1([F[0] and F[k]/F[0] for k in range(N)])
    lF[0] = mp.log(F[0])
    oop = [mp.mpf((-1)**k) for k in range(N)]
    q = smul(oop, lF); Fp = sexp0([mp.mpf(0)] + [-x for x in q[1:]])
    Fp = [x*mp.exp(-q[0]) for x in Fp]
    L, h = e['L'], e['h']; sig = (1-t)*e['s0'] + t*e['s1']
    P2 = smul(Phi, Phi)
    X = [-2*sig*Phi[k]/L + (sig**2 + h**2)*P2[k]/L**2 for k in range(N)]
    sq = [mp.mpf(1)] + [mp.mpf(0)]*(N-1); pw = list(sq)
    for j in range(1, N):
        pw = smul(pw, X); sq = [sq[k] + mp.binomial(mp.mpf(1)/2, j)*pw[k] for k in range(N)]
    S = [L*x for x in sq]
    return smul(Fp, S)

def P_coeffs_t(bs, ws, lams):
    """sum over all edges of int_0^1 E_k dt, k = 0..3 (so P^{(k)}(1) = k! * value)"""
    tot = [mp.mpf(0)]*N
    for e in edges(bs, ws, lams):
        for k in range(N):
            tot[k] += mp.quad(lambda t: E_series(e, t)[k], [0, mp.mpf(1)/2, 1])
    return tot

def zeta3_edge_law(bs, ws):
    """closed-form zeta(3)-coefficient of P'''(1): sum_e 6 al be (h^2/L)(3(al u0 + be u1) - 2)"""
    Z = mp.mpf(0); parts = []
    for e in edges(bs, ws, [0]*len(bs)):
        u0 = -e['s0']/e['L']; u1 = e['s1']/e['L']
        z = 6*e['alpha']*e['beta']*e['h']**2/e['L']*(3*(e['alpha']*u0 + e['beta']*u1) - 2)
        parts.append(z); Z += z
    return Z, parts

def P_direct(bs, ws, lams, p):
    """independent polar perimeter: rho = S_p^{-1/p}, int sqrt(rho^2 + rho'^2), split at kinks"""
    d = p - 1
    def f(th):
        S = mp.mpf(0); Sp = mp.mpf(0)
        for b, w, lam in zip(bs, ws, lams):
            sn = mp.sin(th - b); s = abs(sn)
            W = w*mp.exp(d*lam)
            S += W*s**p
            Sp += W*p*s**(p-1)*mp.cos(th - b)*mp.sign(sn)
        rho = S**(-1/p); drho = -(1/p)*S**(-1/p - 1)*Sp
        return mp.sqrt(rho**2 + drho**2)
    cuts = sorted(set([mp.mpf(0)] + [b % mp.pi for b in bs] + [mp.pi]))
    return 2*mp.quad(f, cuts)

def fd_derivs(bs, ws, lams, h=mp.mpf('0.004'), K=7):
    xs = [k*h for k in range(-K, K+1)]
    ys = [P_direct(bs, ws, lams, 1 + x) for x in xs]
    M = mp.matrix([[x**j for j in range(2*K+1)] for x in xs])
    cf = mp.lu_solve(M, mp.matrix(ys))
    return [cf[0], cf[1], 2*cf[2], 6*cf[3]]
