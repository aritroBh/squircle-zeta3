import mpmath as mp
mp.mp.dps = 60
VS = [(1,1),(1,-1),(-1,1),(-1,-1)]
def hF1(th, tau):
    c, s = mp.cos(th), mp.sin(th); m_ = abs(c)+abs(s)
    L = mp.fsum(mp.exp((v[0]*c+v[1]*s-m_)/tau) for v in VS)
    return m_*mp.exp(tau*mp.log(L))
def tent(th):
    c, s = mp.cos(th), mp.sin(th)
    return abs(c)+abs(s)
def nev(pts):
    n = len(pts); Q = [[mp.mpf(0)]*n for _ in range(n)]
    for i in range(n): Q[i][0] = pts[i][1]
    for j in range(1,n):
        for i in range(n-j):
            xi, xj = pts[i][0], pts[i+j][0]
            Q[i][j] = (xj*Q[i][j-1] - xi*Q[i+1][j-1])/(xj-xi)
    return Q[0][n-1]
for s in [mp.mpf('1'), mp.mpf('2'), mp.mpf('3')]:
    taus = [mp.mpf('0.004')/2**j for j in range(5)]
    vals = []
    for t in taus:
        vals.append(mp.log(hF1(t*s, t)/tent(t*s))/t)
    D1 = nev(list(zip(taus,vals)))
    D2 = nev([(t,(v-D1)/t) for t,v in zip(taus,vals)])
    D3 = nev([(t,(v-D1-D2*t)/t**2) for t,v in zip(taus,vals)])
    g3_pred = (abs(s)**3/6)*(1-mp.tanh(s))
    g1_pred = mp.log(1+mp.exp(-2*abs(s)))
    print(f"s={s}: D1={mp.nstr(D1,16)} (pred {mp.nstr(g1_pred,16)}) D2={mp.nstr(D2,12)} D3={mp.nstr(D3,16)} (pred {mp.nstr(g3_pred,16)})")
