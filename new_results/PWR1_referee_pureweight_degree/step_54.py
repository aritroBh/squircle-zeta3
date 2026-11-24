import mpmath as mp, json
mp.mp.dps = 100
data = json.load(open('quad_rhombus.json'))
taus = [mp.mpf(t) for t in data['taus']]
def fit_ls(ts, vs, kmin=2, kmax=6):
    A = mp.matrix(len(ts), kmax-kmin+1); b = mp.matrix(len(ts), 1)
    for i, t in enumerate(ts):
        for j, k in enumerate(range(kmin, kmax+1)): A[i, j] = t**k
        b[i] = vs[j*0+i]
    return mp.lu_solve(A.T*A, A.T*b)
print(f"{'(a,b)':6} {'w2/((1-q2)z2/4)':18} {'w3/((1-q2)z3/4)':18}  w4, w5 purity below")
for ab, vs_s in data['vals'].items():
    a, b = map(int, ab.split(','))
    q2 = mp.mpf(a*a-b*b)**2/mp.mpf(a*a+b*b)**2
    vs = [mp.mpf(s) for s in vs_s]
    sol = fit_ls(taus, vs)
    w2, w3, w4, w5, w6 = sol
    r2 = w2/((1-q2)*mp.zeta(2)/4); r3 = w3/((1-q2)*mp.zeta(3)/4)
    print(f"({a},{b}):  w2-ratio={mp.nstr(r2,12)}  w3-ratio={mp.nstr(r3,12)}")
    # purity w4: should be rational(q) * pi^4: check residual vs best rational multiple and PSLQ
    x = w4/mp.pi**4
    print(f"    w4/pi^4 = {mp.nstr(x, 18)}  (q^2={mp.nstr(q2,6)})")
    print(f"    w5 = {mp.nstr(w5, 18)}")
json.dump({}, open('tmp.json','w'))
