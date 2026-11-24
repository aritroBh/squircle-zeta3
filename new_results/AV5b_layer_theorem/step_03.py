
import mpmath as mp
mp.mp.dps = 50

def S_p(phi, m, p):
    return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
def rho_p(phi, m, p):
    return S_p(phi, m, p)**(-mp.mpf(1)/p)
def h_p(theta, m, p):
    f = lambda phi: rho_p(phi, m, p)*mp.cos(theta - phi)
    N = 1440
    phis = [2*mp.pi*k/N for k in range(N)]
    vals = [f(ph) for ph in phis]
    k0 = max(range(N), key=lambda k: vals[k])
    try:
        root = mp.findroot(lambda ph: mp.diff(f, ph), phis[k0], tol=mp.mpf('1e-42'))
        return f(root)
    except Exception:
        return vals[k0]

# fan half-width: distance from corner theta0 to adjacent tent corners
for m in [2,3,4]:
    th0 = mp.pi/2 - mp.pi/(2*m)
    # adjacent corners of the tent are at edge normals of adjacent edges: spaced pi/m
    W = mp.pi/(2*m)*mp.mpf('0.6')
    print(f"=== m={m}: fixed-u deviation D(u)=log(h_p/h_1)/tau ===")
    for u in [W*mp.mpf(x) for x in ['0.2','0.5','0.9']]:
        row = []
        for d in ['0.1','0.05','0.025']:
            d = mp.mpf(d); p = 1+d; tau = 1-1/p
            D = mp.log(h_p(th0+u, m, p)/h_p(th0+u, m, 1))/tau
            row.append(mp.nstr(D, 10))
        print(f"  u={mp.nstr(u,6)}: D(tau=.0909,.0476,.0244) = {row}")
    print(f"=== m={m}: scaled profile D(s) near corner ===")
    d = mp.mpf('0.02'); p = 1+d; tau = 1-1/p
    for s in ['0','0.5','1','2','4','8']:
        s = mp.mpf(s); u = tau*s
        D = mp.log(h_p(th0+u, m, p)/h_p(th0+u, m, 1))/tau
        print(f"  s={mp.nstr(s,4)}: D={mp.nstr(D,10)}")
