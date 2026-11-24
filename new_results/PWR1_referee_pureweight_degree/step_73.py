import mpmath as mp, json, pickle
mp.mp.dps = 90
Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m*mp.log1p(mp.e**(-s))**j, [0, 0.5, 2, 8, 32, mp.inf])
    return Icache[(m,j)]
e = {}
with open('e_rm_mine.pkl','rb') as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(',')); e[(r,m)] = mp.mpf(val)
def Dtilde(n, r):
    return 2*mp.fsum(e[(r,n-1-j)]*I(n-1-j,j)/mp.factorial(j) for j in range(1, n-2*r))
def Theta(n, u):
    return mp.fsum(Dtilde(n, m)*(-4*u)**m/mp.factorial(2*m) for m in range(0, n//2))

d = json.load(open('quad_ak.json'))
taus = [mp.mpf(t) for t in d['taus']]
ak = {int(k): [mp.mpf(x) for x in row] for k, row in d['data'].items()}
Ak = {k: -(-1)**k*4*mp.sqrt(2)/(mp.pi*(16*k*k-1)) for k in range(0, 7)}
KMAX = 7  # extract tau^2..tau^8
coef = {}
for k in range(0, 7):
    vs = [v - Ak[k] for v in ak[k]]   # FIX: subtract constant A_k
    A = mp.matrix(len(taus), KMAX); b = mp.matrix(len(taus), 1)
    for i, t in enumerate(taus):
        for j, n in enumerate(range(2, 2+KMAX)): A[i, j] = t**n
        b[i] = vs[i]
    sol = mp.lu_solve(A.T*A, A.T*b)
    coef[k] = {n: (-1)**k*mp.pi/mp.sqrt(2)*sol[n-2] for n in range(2, 2+KMAX)}
print("extracted Theta_n(k^2) vs my D~-built Theta_n(u) (rel. err), and degree test (7 values u=0,1,4,9,16,25,36)")
for n in range(2, 9):
    deg = n//2 - 1
    us = [k*k for k in range(0, 7)]
    vals = [coef[k][n] for k in range(0, 7)]
    rels = [abs(vals[i]-Theta(n, us[i]))/max(abs(Theta(n, us[i])), mp.mpf('1e-30')) for i in range(7)]
    def polyfit(us, vs, dg):
        A = mp.matrix(len(us), dg+1); b = mp.matrix(len(us), 1)
        for i, u in enumerate(us):
            for j in range(dg+1): A[i, j] = mp.mpf(u)**j
            b[i] = vs[i]
        return mp.lu_solve(A.T*A, A.T*b)
    fit = polyfit(us, vals, deg)
    resid = max(abs(vals[i] - mp.fsum(fit[j]*us[i]**j for j in range(deg+1))) for i in range(7))
    over = polyfit(us, vals, deg+1)
    scale = max(abs(v) for v in vals)
    print(f"n={n} (deg {deg}): maxrelerr vs Theta = {mp.nstr(max(rels),2)}; deg-fit max|resid| = {mp.nstr(resid,2)} "
          f"(scale {mp.nstr(scale,2)}); overfit u^{deg+1}/top = {mp.nstr(abs(over[deg+1]*36**(deg+1))/scale,2)}")
json.dump({str(k): {str(n): str(coef[k][n]) for n in coef[k]} for k in coef}, open('theta_extracted.json','w'))
