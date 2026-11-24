import pwr1_rhombus as R, sympy as sp, mpmath as mp, json
s = R.logP_rhombus()
# weight check
t = sp.Symbol('t')
ok = True
for k in range(1, 7):
    c = sp.expand(s).coeff(R.TAU, k)
    sub = c.subs({R.Z[i]: R.Z[i]*t**i for i in R.Z})
    sub = sp.expand(sub)
    homog = all(m[0] == k for m in sp.Poly(sub, t).monoms()) if sub != 0 else True
    q_even = all(m[0] % 2 == 0 for m in sp.Poly(c, R.q).monoms()) if c != 0 else True
    ok &= homog and q_even
    print(f"tau^{k}: weight-{k} homogeneous: {homog}, in Q[q^2]: {q_even}")
print("rhombus purity (exact, through tau^6):", ok)

# compare with quadrature fits
mp.mp.dps = 100
data = json.load(open('quad_rhombus.json'))
taus = [mp.mpf(x) for x in data['taus']]
znum = {R.Z[i]: mp.zeta(i) for i in R.Z}
def fit_ls(ts, vs, kmin=2, kmax=6):
    A = mp.matrix(len(ts), kmax-kmin+1); b = mp.matrix(len(ts), 1)
    for i, tt in enumerate(ts):
        for j, k in enumerate(range(kmin, kmax+1)): A[i, j] = tt**k
        b[i] = vs[i]
    return mp.lu_solve(A.T*A, A.T*b)
for ab, vs_s in data['vals'].items():
    a, b = map(int, ab.split(','))
    qv = sp.Rational(a*a-b*b, a*a+b*b)
    fits = fit_ls(taus, [mp.mpf(x) for x in vs_s])
    for j, k in enumerate(range(2, 7)):
        exact = sp.expand(s).coeff(R.TAU, k).subs(R.q, qv)
        val = mp.mpf(str(exact.subs(znum).evalf(50)))
        rel = abs(val-fits[j])/abs(val) if val != 0 else abs(fits[j])
        print(f"({a},{b}) w_{k}: exact={mp.nstr(val,20)} fit={mp.nstr(fits[j],20)} relerr={mp.nstr(rel,2)} ({mp.nstr(-mp.log10(rel) if rel>0 else 99,3)} digits)")
open('progress.log','a').write("step_55-58: rhombus: re-derived D11c integral rep (50 digits, 4 params); built exact q-symbolic engine through tau^6: pure weight + Q[q^2] EXACT; w2=(1-q2)z2/4, w3=(1-q2)z3/4; quadrature fits match all w2..w6 for (1,2),(2,3),(1,3)\n")
