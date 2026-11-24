"""verify_extend.py <dps> <n> k [k ...]
For k beyond the author's range (k = 11..15; exact pieces re-derived by ME with the author's fastred from WT2 binteg_<k>.pkl, in $VERIFY_EXTRA):
 (i) sum over all weights of my GL integrals vs M_k integrated directly from the squircle definition,
 (ii) W_4(k) (my GL) vs the quadratic c_0 + c_1 k + c_2 k^2 with the c_j from my 800-digit quadrature (verify_const.py halves).
PRECISION: mp.mp.dps set after all imports and printed with every residual."""
import sys, os, json
sys.dont_write_bytecode = True
from verify_lib import *
dps, n = int(sys.argv[1]), int(sys.argv[2]); ks = [int(x) for x in sys.argv[3:]]
set_dps(dps)
print("mp.mp.dps =", mp.mp.dps, flush=True)
OUT = os.environ["VERIFY_OUT"]
def c(w): return sum(mp.mpf(json.load(open(os.path.join(OUT, f"const_{w}_h{h}_dps800.json")))["value"]) for h in (0, 1))
cj = [c("c0"), c("c1"), c("c2")]
X = default_X(dps)
for k in ks:
    pcs = pieces_mp(red_exact(k))
    degs = sorted({i + j for (i, j) in pcs})
    fs = [piece_integrand(pcs, [key for key in pcs if sum(key) == d]) for d in degs]
    vals = integrate(fs, X, n)
    direct = integrate([mk_direct_integrand(k)], X, n, halves=(0,))[0]
    W = dict(zip(degs, vals))
    pred = cj[0] + cj[1] * k + cj[2] * k * k
    print(f"k={k}: M_k direct = {mp.nstr(direct, 20)};  sum of all weights - direct = {mp.nstr(sum(vals) - direct, 4)};  "
          f"W_4(k) = {mp.nstr(W[3], 30)};  W_4(k) - (c0+c1 k+c2 k^2) = {mp.nstr(W[3] - pred, 4)}   (mp.dps {mp.mp.dps})", flush=True)
