"""Referee's own polar jet integrator (no project code imported).
K_p = { sum_j w_j |<u_j,x>|^p <= 1 }  (power convention: w_j = 1, u_j = v_j).
P(p) = int_0^{2pi} r sqrt(1+(r_theta/r)^2) dtheta, r = F(u)^{-1/p}.
The integrand is expanded at each node as a truncated power series in d = p-1
(differentiation under the integral sign), and integrated by tanh-sinh,
split at the kink directions; endpoint terms are evaluated from the endpoint
offset so no cancellation occurs near a corner."""
import mpmath as mp

def s_mul(a, b, N):
    return [mp.fsum(a[k]*b[n-k] for k in range(n+1)) for n in range(N)]

def s_exp(f, N):
    g = [mp.exp(f[0])] + [mp.mpf(0)]*(N-1)
    for n in range(1, N):
        g[n] = mp.fsum(k*f[k]*g[n-k] for k in range(1, n+1))/n
    return g

def s_log(f, N):
    g = [mp.log(f[0])] + [mp.mpf(0)]*(N-1)
    for n in range(1, N):
        g[n] = (f[n] - mp.fsum(k*g[k]*f[n-k] for k in range(1, n))/n)/f[0]
    return g

def s_div(a, b, N):
    q = [mp.mpf(0)]*N
    for n in range(N):
        q[n] = (a[n] - mp.fsum(q[k]*b[n-k] for k in range(n)))/b[0]
    return q

def s_sqrt(f, N):
    g = [mp.sqrt(f[0])] + [mp.mpf(0)]*(N-1)
    for n in range(1, N):
        g[n] = (f[n] - mp.fsum(g[k]*g[n-k] for k in range(1, n)))/(2*g[0])
    return g

def integrand(theta, da, db, ja, jb, U, W, N):
    c, s = mp.cos(theta), mp.sin(theta)
    F = [mp.mpf(0)]*N; G = [mp.mpf(0)]*N
    for j, (ux, uy) in enumerate(U):
        if j == ja:
            nr = mp.sqrt(ux*ux+uy*uy); a = nr*mp.sin(da); ap = nr*mp.cos(da)
        elif j == jb:
            nr = mp.sqrt(ux*ux+uy*uy); a = nr*mp.sin(db); ap = -nr*mp.cos(db)
        else:
            v = ux*c+uy*s; vp = -ux*s+uy*c
            a = abs(v); ap = vp if v > 0 else -vp
        la = mp.log(a); term = mp.mpf(1)
        for k in range(N):
            F[k] += W[j]*a*term; G[k] += W[j]*ap*term
            term = term*la/(k+1)
    inv1pd = [(-1)**k for k in range(N)]
    logF = s_log(F, N)
    r = s_exp([-x for x in s_mul(logF, inv1pd, N)], N)
    Q = s_div(G, F, N)
    Q2 = s_mul(Q, Q, N); Q2[0] += 1
    return s_mul(r, s_sqrt(Q2, N), N)

def kinks(U):
    out = []
    for j, (ux, uy) in enumerate(U):
        phi = mp.atan2(ux, -uy)
        phi = phi % mp.pi
        out.append((phi, j))
    out.sort()
    return out

def perimeter_jet(U, W=None, N=4, h=mp.mpf(1)/64, cut=None):
    """returns list [P_0..P_{N-1}] with P(1+d) = sum P_k d^k (asymptotic)."""
    if W is None: W = [mp.mpf(1)]*len(U)
    if cut is None: cut = mp.mpf(10)**(-(mp.mp.dps+5))
    ks = kinks(U)
    m = len(ks)
    tot = [mp.mpf(0)]*N
    for i in range(m):
        tha, ja = ks[i]
        thb, jb = ks[(i+1) % m]
        if i == m-1: thb = thb + mp.pi
        Wd = thb - tha
        acc = [mp.mpf(0)]*N
        k = 0
        while True:
            done = True
            for t in ([k*h] if k == 0 else [k*h, -k*h]):
                q = mp.exp(-mp.pi*mp.sinh(t))
                s = 1/(1+q); sc = q/(1+q)
                if min(s, sc) < cut: continue
                done = False
                wgt = mp.pi*mp.cosh(t)*s*sc*Wd*h
                da = Wd*s; db = Wd*sc
                th = tha+da if s < sc else thb-db
                f = integrand(th, da, db, ja, jb, U, W, N)
                for n in range(N): acc[n] += wgt*f[n]
            if done and k > 0: break
            k += 1
        for n in range(N): tot[n] += 2*acc[n]
    return tot

def perimeter_direct(U, p, W=None, h=mp.mpf(1)/64):
    """direct P(p) (no expansion), same quadrature."""
    if W is None: W = [mp.mpf(1)]*len(U)
    p = mp.mpf(p)
    ks = kinks(U); m = len(ks); tot = mp.mpf(0)
    cut = mp.mpf(10)**(-(mp.mp.dps+5))
    for i in range(m):
        tha, ja = ks[i]; thb, jb = ks[(i+1) % m]
        if i == m-1: thb += mp.pi
        Wd = thb-tha; acc = mp.mpf(0); k = 0
        while True:
            done = True
            for t in ([k*h] if k == 0 else [k*h, -k*h]):
                q = mp.exp(-mp.pi*mp.sinh(t)); s = 1/(1+q); sc = q/(1+q)
                if min(s, sc) < cut: continue
                done = False
                wgt = mp.pi*mp.cosh(t)*s*sc*Wd*h
                da = Wd*s; db = Wd*sc
                th = tha+da if s < sc else thb-db
                c, sn = mp.cos(th), mp.sin(th)
                F = 0; G = 0
                for j, (ux, uy) in enumerate(U):
                    if j == ja:
                        nr = mp.sqrt(ux*ux+uy*uy); a = nr*mp.sin(da); ap = nr*mp.cos(da)
                    elif j == jb:
                        nr = mp.sqrt(ux*ux+uy*uy); a = nr*mp.sin(db); ap = -nr*mp.cos(db)
                    else:
                        v = ux*c+uy*sn; vp = -ux*sn+uy*c; a = abs(v); ap = vp if v > 0 else -vp
                    F += W[j]*a**p; G += W[j]*p*a**(p-1)*ap
                r = F**(-1/p); Q = G/(p*F)
                acc += wgt*r*mp.sqrt(1+Q*Q)
            if done and k > 0: break
            k += 1
        tot += 2*acc
    return tot

def vertices(U, W=None):
    """vertices of K_1 in ccw order (one per kink direction, both signs)."""
    if W is None: W = [mp.mpf(1)]*len(U)
    out = []
    for phi, j in kinks(U):
        for ph in (phi, phi+mp.pi):
            d = (mp.cos(ph), mp.sin(ph))
            nrm = mp.fsum(W[i]*abs(U[i][0]*d[0]+U[i][1]*d[1]) for i in range(len(U)) if i != j)
            out.append((ph, (d[0]/nrm, d[1]/nrm)))
    out.sort(key=lambda z: z[0])
    return [v for _, v in out]
