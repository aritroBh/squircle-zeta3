
code = r'''
import json, os
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, cos, sin, exp, log, nstr, gamma

mp.dps = 120
raw = json.load(open("../D11b2_stokes_jump/step2_num450_N180.json"))
f = [mpf(v) for v in raw["f450"]]
Nmax = 180
assert len(f) >= Nmax+1

# Borel coefficients b_n = f_n/n!
b = [f[n]/factorial(n) for n in range(Nmax+1)]

A = mpc(0, pi/2)          # singularity at +i pi/2
invA = 2/(mpc(0,1)*pi)    # A^{-1} = 2/(i pi)

# ---- 1) exponent gamma from period-4 ratios: b_{n+4}/b_n = (2/pi)^4 ((n+4)/n)^{-(gamma+1)} (1+O(1/n^2))
#     e_n := log( (pi/2)^4 b_{n+4}/b_n ) / log( n/(n+4) )  -> gamma+1 = 3/2
print("exponent estimates e_n -> gamma+1 (expect 1.5 if branch (1-t/A)^{1/2}):")
es = {}
for n in range(20, Nmax-4):
    r = (pi/2)**4 * b[n+4]/b[n]
    es[n] = log(r)/log(mpf(n)/(n+4))
for n in [40, 80, 120, 160, 172]:
    print("  n=%3d: e_n = %s" % (n, nstr(es[n], 12)))
# Richardson in 1/n: e_n = g1 + c1/n + c2/n^2 + ...  fit on windows
def richardson(seq, ns, order):
    # fit seq[n] = s + sum_{j=1..order} c_j n^{-j} by least squares
    m = order+1
    M = mp.matrix(len(ns), m); y = mp.matrix(len(ns),1)
    for i,n in enumerate(ns):
        for j in range(m):
            M[i,j] = mpf(n)**(-j)
        y[i] = seq[n]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G, rhs)
    return sol[0], sqrt(abs((M*sol-y).T*(M*sol-y))[0]/len(ns))
print("\nRichardson extrapolants for gamma+1:")
for lo, hi, order in [(100,172,2),(120,172,3),(140,172,3),(150,172,4)]:
    ns = list(range(lo, hi))
    s, rms = richardson(es, ns, order)
    print("  window [%d,%d] order %d: gamma+1 = %s  (rms %s)" % (lo,hi,order,nstr(s,14), nstr(rms,2)))

# ---- 2) amplitude: d_n = b_n (pi/2)^n n^{3/2} / cos-pattern should -> sqrt2 (with 1/n corrections)
print("\namplitude d_n = b_n (pi/2)^n n^{3/2} / (sqrt2 cos(pi(3-2n)/4)) -> 1:")
pat = lambda n: cos(pi*(3-2*n)/4)
ds = {}
for n in range(6, Nmax+1):
    if abs(pat(n)) < 0.5: continue
    ds[n] = b[n]*(pi/2)**n * n**mpf('1.5') / (sqrt(2)*pat(n))
for n in [40, 80, 120, 160, 180]:
    print("  n=%3d: d_n = %s" % (n, nstr(ds[n], 12)))
print("\nRichardson extrapolants for amplitude (expect 1):")
for lo, hi, order in [(100,180,2),(120,180,3),(140,180,4)]:
    ns = list(range(lo, hi))
    s, rms = richardson(ds, ns, order)
    print("  window [%d,%d] order %d: amp = %s  (rms %s)" % (lo,hi,order,nstr(s,14), nstr(rms,2)))

# ---- 3) implied Borel branch constant c_+ in B(t) ~ c_+ (1 - t/A)^{1/2} and Stokes constant
# b_n ~ c_+ A^{-n} [z^n](1-z)^{1/2} + c.c ;  [z^n](1-z)^{1/2} ~ -n^{-3/2}/(2 sqrt pi)
# fit: b_n n^{3/2} = -(1/(2 sqrt pi))(c_+ A^{-n} + c_- A^{-n}) ; c_- = conj(c_+)
# with c_+ = -sqrt(2 pi) e^{3 i pi/4}: prediction (from D11): check
cplus_pred = -sqrt(2*pi)*exp(mpc(0, 3*pi/4))
print("\npredicted c_+ = -sqrt(2 pi) e^{3 i pi/4} =", nstr(cplus_pred, 15))
# numerical c_+ from d_n: sqrt2 cos = (1/sqrt2)(e^{3i pi/4}(2/(i pi))^n + e^{-3i pi/4}(-2/(i pi))^n)
# => c_+ A^{-n}(-1/(2 sqrt pi)) = (1/sqrt2) e^{3i pi/4} (2/(i pi))^n  => c_+ = -2 sqrt pi /sqrt2 e^{3i pi/4} = -sqrt(2 pi) e^{3i pi/4}
cpl = -2*sqrt(pi)/sqrt(2)*exp(mpc(0,3*pi/4))
print("re-derived c_+ =", nstr(cpl, 15))
# Stokes/trans-monomial constant: S such that F ⊃ S eta^{1/2} e^{-A/eta}: S = c_+ * 2 i Gamma(3/2) * phase
S_mag = 2*sqrt(pi)
print("|S| = 2 sqrt pi =", nstr(S_mag, 15), " ; D11 jump constant C = 2 sqrt pi =", nstr(2*sqrt(pi),15))
# local constant gamma_A in B_F(A+zeta) ~ gamma_A zeta^{1/2}: gamma_A = c_+ (-1/A)^{1/2}
gA = cpl*sqrt(-1/A)
print("gamma_A (A=+i pi/2): predicted 2 : computed", nstr(gA, 12))
'''
open("step_04.py","w").write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_04.py"],capture_output=True,text=True)
print(r.stdout); print("STDERR:", r.stderr[-1500:] if r.stderr else "")
