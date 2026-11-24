import os
os.chdir("<project>/research/08_theory/AV6_squircle_lattice_zeta")
exec(open("squircle.py").read().replace("def setdps(d): mp.dps = d","def setdps(d): mp.dps = d"))
setdps(60)
from mpmath import mp, mpf, nstr, log, power, fsum, zeta
import mpmath
s = 4
# G_x(4) = sum_k k^{-5} sigma_x(k), sigma_x(k) = sum_{m<k} m ln m (k-m) ln(k-m)
K = 600
P = mpf(0)
for k in range(2, K+1):
    sx = fsum([m*log(m)*(k-m)*log(k-m) for m in range(1, k)])
    P += power(k, -5)*sx
# tail: sigma_x(k) ~ k^3(ln^2 k/6 - 5 ln k/18 + 1/18) + k^2 ln^2 k/2  (EM next term)
# residual after lead: fit from exact shells at K-30..K as (a k^2 ln^2 k + b k^2 ln k)
xs, ys = [], []
for k in range(K-30, K+1):
    sx = fsum([m*log(m)*(k-m)*log(k-m) for m in range(1, k)])
    lead = k**3*(log(k)**2/6 - 5*log(k)/18 + mpf(1)/18)
    xs.append([k**2*log(k)**2, k**2*log(k)]); ys.append(sx-lead)
A = mpmath.matrix(len(xs), 2); b = mpmath.matrix(len(xs), 1)
for i in range(len(xs)):
    A[i,0] = xs[i][0]; A[i,1] = xs[i][1]; b[i,0] = ys[i]
c = mpmath.lu_solve(A.T*A, A.T*b)
tail = mpf(0)
for k in range(K+1, 40*K+1):
    model = c[0,0]*k**2*log(k)**2 + c[1,0]*k**2*log(k)
    lead = k**3*(log(k)**2/6 - 5*log(k)/18 + mpf(1)/18)
    tail += power(k, -5)*(lead + model)
print("G_x(4) improved:", nstr(P+tail, 25))
print("needed for identity: ~0.13397 ; code had: 0.13094")
# also G_x(5), G_x(6) quickly for the c2 fits at s=5,6? skip.
