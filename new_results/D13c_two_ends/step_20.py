
code = r'''
import json
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, cos, log, nstr, exp
mp.dps = 120
raw = json.load(open("../D11b2_stokes_jump/step2_num450_N180.json"))
f = [mpf(v) for v in raw["f450"]]
Nmax = 180
b = [f[n]/factorial(n) for n in range(Nmax+1)]
pat = lambda n: cos(pi*(3-2*n)/4)

# ---- exponent: e_n = gamma+1 + c/n fit on windows (2-param exact LSQ via normal equations, well conditioned)
es = {}
for n in range(20, Nmax-3):
    r = (pi/2)**4 * b[n+4]/b[n]
    es[n] = log(r)/log(mpf(n)/(n+4))
def fit2(xs, ys):
    n = len(xs); sx = sum(xs); sxx = sum(x*x for x in xs); sy = sum(ys); sxy = sum(x*y for x,y in zip(xs,ys))
    d = n*sxx-sx*sx
    return (sxx*sy-sx*sxy)/d, (n*sxy-sx*sy)/d
print("fit e_n = g + c/n over windows:")
for lo,hi in [(100,172),(120,172),(140,172),(150,172),(160,172)]:
    xs = [1/mpf(n) for n in range(lo,hi+1)]; ys = [es[n] for n in range(lo,hi+1)]
    g,c = fit2(xs,ys)
    print("  [%3d,%3d]: g = %s  (|g-1.5|=%s)" % (lo,hi,nstr(g,13),nstr(abs(g-mpf('1.5')),3)))

# ---- amplitude+corrections: r_n = b_n / [sqrt2 (2/pi)^n cos n^{-3/2}] = 1 + c1/n + c2/n^2 + ...
r_ = {}
for n in range(6,Nmax+1):
    r_[n] = b[n]/(sqrt(2)*(2/pi)**n*pat(n)*mpf(n)**mpf('-1.5'))
def lsq(ns, deg):
    m = deg+1
    M = mp.matrix(len(ns), m); y = mp.matrix(len(ns),1)
    for i,n in enumerate(ns):
        for j in range(m): M[i,j] = mpf(n)**(-j)
        y[i] = r_[n]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs)
    res = M*sol-y
    rms = sqrt(sum(res[i]**2 for i in range(len(ns)))/len(ns))
    return sol, rms
print("\nleast-squares r_n = 1 + c1/n + ... : constant and c1 (analytic c1 prediction: see below)")
for lo,deg in [(60,3),(80,4),(100,4),(100,5),(120,5)]:
    ns = list(range(lo,181))
    sol,rms = lsq(ns,deg)
    print("  [%3d,180] deg %d: const=%s  c1=%s  rms=%s" % (lo,deg,nstr(sol[0],14),nstr(sol[1],8),nstr(rms,3)))

# analytic c1: Borel branch B_F(A+zeta) ~ 2 zeta^{1/2}(1 + d1 zeta + ...) with d1 from
# g(x) = e^{-x/2}(cosh x)^{1/2} near x = i pi/2:  g ~ zeta^{1/2}(1 - zeta/2 + ...) [from cosh(i pi/2+zeta)=i sinh zeta, e^{-x/2}]
# coefficient ratio d1 = -1/2 -> b_n correction: (1 - t/A)^1 extra power shifts exponent: [z^n](1-z)^{3/2}/[z^n](1-z)^{1/2} ~ (3/2)/n ... 
# predict c1/n coefficient: [z^n](1-z)^{1/2}(1+d1'(1-z)) : ratio term = d1' * [z^n](1-z)^{3/2}/[z^n](1-z)^{1/2} = d1' * (3/2)/(n+1/2) * ...
# more precisely: [z^n](1-z)^{3/2} ~ -(3/(4 sqrt pi)) n^{-5/2}; [z^n](1-z)^{1/2} ~ -(1/(2 sqrt pi)) n^{-3/2}: ratio = (3/2)/n
# so c1 = d1' * 3/2 with d1' = d1*(A factor)... compute numerically below from the fit.
'''
open("step_06.py","w").write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_06.py"],capture_output=True,text=True)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.stderr else "")
