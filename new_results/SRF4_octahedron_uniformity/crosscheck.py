"""Cross-checks of the coefficient engine: paper thm:oct log-series tau^2..tau^6, SRF3 K120 data (all digits),
precision stability K150 (dps 330) vs K80 (dps 160/220)."""
import json
from mpmath import mp, mpf, pi, zeta, log, nstr, factorial
mp.dps = 300
A = json.load(open('coeffs_K150_dps330.json')); a = [mpf(x) for x in A['a_tau']]; K = A['K']
# log((2-p)S/(4sqrt3)) in tau: 2-p = (1-2tau)/(1-tau)
def mul(x, y, n): return [sum(x[i]*y[k-i] for i in range(k+1)) for k in range(n)]
n = 8; geo = [mpf(1)]*n                       # 1/(1-tau)
f = mul(mul([mpf(1), mpf(-2)]+[mpf(0)]*(n-2), geo, n), a[:n], n)
L = [mpf(0)]*n; L[0] = log(f[0])
for k in range(1, n): L[k] = (k*f[k]-sum(j*L[j]*f[k-j] for j in range(1, k)))/(k*f[0])
ref = {2: -pi**2/9, 3: -16*zeta(3)/9, 4: -37*pi**4/810, 5: -160*zeta(5)/27, 6: 16*zeta(3)**2/81-572*pi**6/76545}
print("log-series vs thm:oct:", {k: nstr(L[k]-ref[k], 3) for k in ref}, " L1 =", nstr(L[1], 3))
s = json.load(open('../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json'))
dig = min(-log(abs(mpf(s['eta'][k])-a[k])/abs(a[k]), 10) if mpf(s['eta'][k]) != a[k] else 999 for k in range(121))
digd = min(-log(abs(mpf(s['d'][k])-mpf(A['b_d'][k]))/abs(mpf(A['b_d'][k])), 10) if mpf(s['d'][k]) != mpf(A['b_d'][k]) else 999 for k in range(121))
print("agreement with SRF3 K120 file (tau, d), min relative digits over k<=120:", nstr(dig, 5), nstr(digd, 5))
B = json.load(open('coeffs_K80_dps220.json'))
print("K150/dps330 vs K80/dps220, min relative digits k<=80:", nstr(min(-log(abs(mpf(B['a_tau'][k])-a[k])/abs(a[k]), 10) if mpf(B['a_tau'][k]) != a[k] else 999 for k in range(81)), 5))
print("a_k, k=0..6:", [nstr(a[k], 12) for k in range(7)])
