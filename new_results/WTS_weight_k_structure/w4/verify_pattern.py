"""verify_pattern.py <dps> <n> <kmax>
Tests the author's side claim "weight w obeys the polynomial law (degree w-2 in k) for k >= w and fails at k = w-1" with MY quadrature
(verify_lib) on the exact pieces red_k.pkl (k = 1..kmax; k=1 piece file made by the author's w4_k1_probe.py).
Prints W_w(k) for w = 2..6 and the (w-1)-th finite differences Delta^{w-1} W_w(k0) for each window start k0.
PRECISION: mp.mp.dps set after all imports and printed."""
import sys, os, time
sys.dont_write_bytecode = True
from verify_lib import *

dps, n, kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
set_dps(dps)
print("mp.mp.dps =", mp.mp.dps, flush=True)
X = default_X(dps)
WMAX = 6
W = {}                                   # W[(w,k)]
for k in range(1, kmax + 1):
    pcs = pieces_mp(red_exact(k))
    degs = [d for d in range(0, WMAX) if any(i + j == d for (i, j) in pcs)]
    fs = [piece_integrand(pcs, [key for key in pcs if sum(key) == d]) for d in degs]
    vals = integrate(fs, X, n) if fs else []
    for d, v in zip(degs, vals):
        W[(d + 1, k)] = v
    print(f"k={k}: weights present {[d+1 for d in degs]}  (max weight in file {max(i+j for (i,j) in pcs)+1})", flush=True)
for w in range(2, WMAX + 1):
    print(f"--- weight {w}: law = polynomial of degree {w-2} in k; Delta^{w-1} over window k0..k0+{w-1}")
    ks = [k for k in range(1, kmax + 1) if (w, k) in W]
    print("   W_%d(k) = %s" % (w, [(k, mp.nstr(W[(w, k)], 12)) for k in ks]))
    for k0 in range(1, kmax - (w - 1) + 1):
        win = [W.get((w, k0 + t), mp.mpf(0)) for t in range(w)]
        # (w-1)-th forward difference
        d = sum((-1) ** (w - 1 - t) * mp.binomial(w - 1, t) * win[t] for t in range(w))
        print(f"   k0={k0}: Delta^{w-1} W_{w} = {mp.nstr(d, 8)}   (window k={k0}..{k0+w-1};  mp.dps {mp.mp.dps})")
