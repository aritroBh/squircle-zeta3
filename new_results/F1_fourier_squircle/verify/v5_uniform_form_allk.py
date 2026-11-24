"""
VERIFIER V5: Proposition 3b (uniform boundary-layer form) over ALL k = 1..400, not just a sample.
a_k - A_k by composite Gauss-Legendre in theta (double precision, 4000 panels x 24 nodes),
compared with (-1)^k (sqrt2/pi)[1/(4k^2) - pi tau/(2k sinh(2 pi k tau))].
"""
import numpy as np
x, w = np.polynomial.legendre.leggauss(24)
NP = 4000
edges = np.linspace(0, np.pi / 4, NP + 1)
a, b = edges[:-1, None], edges[1:, None]
TH = ((b - a) / 2 * x[None, :] + (a + b) / 2).ravel()
W = ((b - a) / 2 * w[None, :]).ravel()
def amA_all(tau, K):
    q = 1 / tau
    base = np.cos(TH) * np.expm1(tau * np.log1p(np.tan(TH) ** q)) * W
    ks = np.arange(0, K + 1)
    return 8 / np.pi * (np.cos(4 * np.outer(ks, TH)) @ base)
K = 400
for tau in [0.02, 0.01, 0.005]:
    d = amA_all(tau, K)
    k = np.arange(1, K + 1)
    lead = (-1.0) ** k * np.sqrt(2) / np.pi * (1 / (4 * k ** 2) - np.pi * tau / (2 * k * np.sinh(2 * np.pi * k * tau)))
    err = d[1:] - lead
    i = np.argmax(np.abs(err))
    print(f"tau={tau}: max_(1<=k<=400) |a_k - A_k - lead| = {abs(err[i]):.4e} at k={k[i]};  /tau^3 = {abs(err[i])/tau**3:.4f};"
          f"  k=0 check: (-1)^0 pi (a0-A0)/(sqrt2 tau^2) = {np.pi*d[0]/(np.sqrt(2)*tau**2):.10f} (zeta2+zeta3 tau ~ {1.6449340668+1.2020569*tau:.10f})")
