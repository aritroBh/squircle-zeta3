
import json
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, cos, log, nstr, exp
mp.dps = 120
raw = json.load(open("../D11b2_stokes_jump/step2_num450_N180.json"))
f = [mpf(v) for v in raw["f450"]]
Nmax = 180
b = [f[n]/factorial(n) for n in range(Nmax+1)]
pat = lambda n: cos(pi*(3-2*n)/4)
r_ = {}
for n in range(6,Nmax+1):
    r_[n] = b[n]/(sqrt(2)*(2/pi)**n*pat(n)*mpf(n)**mpf('-1.5'))
print("n, n*(r_n-1) for n=170..180 (check period-4 in correction):")
for n in range(170,181):
    print("  n=%d mod4=%d: n*(r_n-1) = %s" % (n, n%4, nstr(n*(r_[n]-1), 10)))

def lsq_on(ns, deg, seq):
    m = deg+1
    M = mp.matrix(len(ns), m); y = mp.matrix(len(ns),1)
    for i,n in enumerate(ns):
        for j in range(m): M[i,j] = mpf(n)**(-j)
        y[i] = seq[n]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs)
    res = M*sol-y
    return sol, sqrt(sum(res[i]**2 for i in range(len(ns)))/len(ns))

print("\nquadrisected amplitude fits r_n = c0 + c1/n + c2/n^2 (per mod-4 class):")
for rho in range(4):
    ns = [n for n in range(100+((rho-100)%4),181,4)]
    sol,rms = lsq_on(ns, 2, r_)
    print("  class %d: c0 = %s  c1 = %s  rms = %s" % (rho, nstr(sol[0],13), nstr(sol[1],8), nstr(rms,2)))
print("\nquadrisected amplitude fits deg 3 [120,180]:")
c1s = {}
for rho in range(4):
    ns = [n for n in range(120+((rho-120)%4),181,4)]
    sol,rms = lsq_on(ns, 3, r_)
    c1s[rho] = sol[1]
    print("  class %d: c0 = %s  c1 = %s  c2 = %s  rms = %s" % (rho, nstr(sol[0],13), nstr(sol[1],8), nstr(sol[2],8), nstr(rms,2)))
# exponent per class too
es = {}
for n in range(20, Nmax-3):
    rr = (pi/2)**4 * b[n+4]/b[n]
    es[n] = log(rr)/log(mpf(n)/(n+4))
print("\nexponent fits e_n = g + c/n + d/n^2 per mod-4 class [100,172]:")
for rho in range(4):
    ns = [n for n in range(100+((rho-100)%4),173,4)]
    sol,rms = lsq_on(ns, 2, es)
    print("  class %d: g = %s  (|g-1.5|=%s)  rms=%s" % (rho, nstr(sol[0],13), nstr(abs(sol[0]-mpf('1.5')),3), nstr(rms,2)))
