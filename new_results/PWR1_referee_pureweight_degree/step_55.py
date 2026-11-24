import mpmath as mp, json
mp.mp.dps = 100
data = json.load(open('quad_rhombus.json'))
taus = [mp.mpf(t) for t in data['taus']]
def fit_ls(ts, vs, kmin=2, kmax=6):
    A = mp.matrix(len(ts), kmax-kmin+1); b = mp.matrix(len(ts), 1)
    for i, t in enumerate(ts):
        for j, k in enumerate(range(kmin, kmax+1)): A[i, j] = t**k
        b[i] = vs[i]
    return mp.lu_solve(A.T*A, A.T*b)
z = mp.zeta
for ab, vs_s in data['vals'].items():
    a, b = map(int, ab.split(','))
    vs = [mp.mpf(s) for s in vs_s]
    w2,w3,w4,w5,w6 = fit_ls(taus, vs)
    r4 = mp.pslq([w4, mp.pi**4, z(2), mp.euler, mp.mpf(1)], tol=mp.mpf('1e-7'), maxcoeff=10**6, maxsteps=10000)
    r5 = mp.pslq([w5, z(5), z(2)*z(3), z(3), z(2), mp.euler, mp.mpf(1)], tol=mp.mpf('1e-6'), maxcoeff=10**6, maxsteps=10000)
    print(f"({a},{b}): pslq w4 [w4,pi^4,z2,gamma,1] = {r4}")
    print(f"       pslq w5 [w5,z5,z2z3,z3,z2,gamma,1] = {r5}")
