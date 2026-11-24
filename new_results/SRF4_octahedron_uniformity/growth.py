"""Compare computed a_k (tau) and b_n (d) with the growth laws of Theorem G / Corollary Gd."""
import json, sys
from mpmath import mp, mpf, mpc, log, pi, sqrt, factorial, re, im, arg, nstr, matrix, lu_solve, cos, sin
mp.dps = 60
d = json.load(open(sys.argv[1])); K = d["K"]
a = [mpf(x) for x in d["a_tau"]]; b = [mpf(x) for x in d["b_d"]]
z0 = mpc(-log(2), pi)/2; R = abs(z0); th = arg(z0)
cT, cD = sqrt(3*pi/8), sqrt(3*pi/4)
lawT = lambda k: -cT*k**mpf(-1.5)*re(z0**(mpf(1)/2-k))
lawD = lambda k: -cD*k**mpf(-1.5)*im(z0**(mpf(1)/2-k))
scale = lambda k: k**mpf(-1.5)*R**(-k)
print("1/R = %s" % nstr(1/R, 16))
print(" k   |a_k/k!|^(1/k)   (a_k/k!-law)/(k^-3/2 R^-k)  k*that   |b_n/n!|^(1/n)  (b_n/n!-lawD)/scale  k*that")
for k in [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120, 150]:
    if k > K: break
    al = a[k]/factorial(k); bl = b[k]/factorial(k)
    r1 = (al-lawT(k))/scale(k); r2 = (bl-lawD(k))/scale(k)
    print("%3d  %s  %s  %s   %s  %s  %s" % (k, nstr(abs(al)**(mpf(1)/k), 10), nstr(r1, 6), nstr(k*r1, 6),
          nstr(abs(bl)**(mpf(1)/k), 10), nstr(r2, 6), nstr(k*r2, 6)))
# limsup proxy: max over windows of |a_k/k!|^(1/k) * R and the k^(-3/2)-corrected root
print("max_{k in [K-20,K]} |a_k/k!|^(1/k)*R =", nstr(max(abs(a[k]/factorial(k))**(mpf(1)/k) for k in range(K-20, K+1))*R, 10))
print("max_{k in [K-20,K]} |a_k/(k! k^-3/2)|^(1/k)*R =", nstr(max(abs(a[k]/factorial(k)/k**mpf(-1.5))**(mpf(1)/k) for k in range(K-20, K+1))*R, 10))
# least-squares fit: a_k/k! / (k^-3/2 R^-k) = A cos(k th) + B sin(k th) + (A1 cos + B1 sin)/k + (A2 cos+B2 sin)/k^2
def fit(seq, k0, k1, ncorr=3):
    rows = []; rhs = []
    for k in range(k0, k1+1):
        y = seq[k]/factorial(k)/scale(k)
        rows.append([f(k*th)/mpf(k)**j for j in range(ncorr) for f in (cos, sin)]); rhs.append(y)
    M = matrix(rows); Y = matrix(rhs)
    sol = lu_solve(M.T*M, M.T*Y); return sol[0], sol[1]
for name, seq, amp_th, ph_th in [("tau", a, sqrt(3*pi*R/8), pi-th/2), ("d", b, sqrt(3*pi*R/4), None)]:
    for (k0, k1) in [(40, K), (70, K)]:
        A, B = fit(seq, k0, k1)
        amp = sqrt(A**2+B**2); ph = mp.atan2(-B, A)   # A cos + B sin = amp cos(k th + ph)
        print("%s-series fit k=%d..%d: amplitude %s (theory %s), phase %s" % (name, k0, k1, nstr(amp, 12), nstr(amp_th, 12), nstr(ph, 12)))
print("theory tau phase pi-theta/2 = %s ; d-series theory: sqrt(3piR/4) sin((k-1/2)theta) = amp cos(k th + phase) with phase = %s" % (
      nstr(pi-th/2, 12), nstr(-th/2-pi/2, 12)))
