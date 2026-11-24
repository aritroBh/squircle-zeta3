import mpmath as mp
mp.mp.dps = 70
VS = [(1,1),(1,-1),(-1,1),(-1,-1)]
def hF1(th, tau):
    c, s = mp.cos(th), mp.sin(th); m_ = abs(c)+abs(s)
    L = mp.fsum(mp.exp((v[0]*c+v[1]*s-m_)/tau) for v in VS)
    return m_*mp.exp(tau*mp.log(L))
def tent(th):
    c, s = mp.cos(th), mp.sin(th)
    return abs(c)+abs(s)
def Q(tau):
    W = 30*tau
    return mp.quad(lambda u: hF1(u, tau) - tent(u), [-W, 0, W])
def nev(pts):
    n = len(pts); Q_ = [[mp.mpf(0)]*n for _ in range(n)]
    for i in range(n): Q_[i][0] = pts[i][1]
    for j in range(1,n):
        for i in range(n-j):
            xi, xj = pts[i][0], pts[i+j][0]
            Q_[i][j] = (xj*Q_[i][j-1] - xi*Q_[i+1][j-1])/(xj-xi)
    errs = [abs(Q_[0][k]-Q_[0][k-1]) for k in range(1,n)]
    best = min(range(1,n), key=lambda k: errs[k-1])
    return Q_[0][best], errs[best-1]
taus = [mp.mpf('0.003')/2**k for k in range(7)]
Qs = [Q(t) for t in taus]
c2, e2 = nev([(t, v/t**2) for t, v in zip(taus, Qs)])
rem = [(t, v - c2*t**2) for t, v in zip(taus, Qs)]
c3, e3 = nev([(t, v/t**3) for t, v in rem])
rem = [(t, v - c3*t**3) for t, v in rem]
c4, e4 = nev([(t, v/t**4) for t, v in rem])
print("c2 =", mp.nstr(c2, 25), " pi2/12 =", mp.nstr(mp.pi**2/12, 25))
print("c3 =", mp.nstr(c3, 25), " zeta(3)/2 =", mp.nstr(mp.zeta(3)/2, 25))
print("c4 =", mp.nstr(c4, 25), "  4*c4 =", mp.nstr(4*c4, 25), " zeta(4)/4 =", mp.nstr(mp.zeta(4)/4, 25))
print("resid 4c4 - zeta(4)/4:", mp.nstr(4*c4 - mp.zeta(4)/4, 6))
