"""verify_quad.py <dps> <n> <wmax> <k> [slow_pickle_dir]
Independent re-evaluation of the weight pieces of M_k with MY Gauss-Legendre quadrature (verify_lib.py) and
(a) the sum over ALL weights compared with M_k integrated DIRECTLY from the squircle definition (no v-substitution, no IBP, no
    reducer at all: [eps^k] of t^{eps-1}+(1-t)^{eps-1}-sqrt(t^{2eps-2}+(1-t)^{2eps-2}) by power-series arithmetic);
(b) W_4(k) compared with the author's stored 120-digit value results/w4_k<k>_dps120.json (read only);
(c) if slow_pickle_dir/slow_red_<k>.pkl exists (original sympy reducer ibp_reduce.reduce_k): pole-by-pole coefficient comparison with the
    fast exact pieces, then the same integrals from the slow pieces.
PRECISION: mp.mp.dps is set after all imports and printed."""
import sys, os, json, time, pickle
sys.dont_write_bytecode = True
from verify_lib import *
import sympy as sp

dps, n, wmax, k = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
slow_dir = sys.argv[5] if len(sys.argv) > 5 else None
set_dps(dps)
print("mp.mp.dps =", mp.mp.dps, " n =", n, " wmax =", wmax, " k =", k, flush=True)
wmax = mp.mpf(wmax)
X = default_X(dps)
t0 = time.time()

red = red_exact(k)
pcs = pieces_mp(red)
degs = sorted({i + j for (i, j) in pcs})
print("fast pieces: keys", sorted(pcs), "degrees", degs)
fs = [piece_integrand(pcs, [key for key in pcs if sum(key) == d]) for d in degs]
fast = integrate(fs, X, n)
Wfast = dict(zip(degs, fast))
for d in degs:
    print(f"  fast  W_{d+1}({k}) = {mp.nstr(Wfast[d], 45)}  (mp.dps {mp.mp.dps})", flush=True)
tot = sum(fast)
print(f"  [{time.time()-t0:.0f}s]")
direct = integrate([mk_direct_integrand(k)], X, n, halves=(0,))[0]
print(f"sum_d W_(d+1)({k}) [fast pieces, my GL]  = {mp.nstr(tot, 45)}")
print(f"M_{k} direct (squircle definition)       = {mp.nstr(direct, 45)}")
print(f"RESIDUAL  sum of all weights - direct M_{k} = {mp.nstr(tot - direct, 5)}   (mp.dps {mp.mp.dps})", flush=True)

fn = os.path.join(W4, "results", f"w4_k{k}_dps120.json")
if os.path.exists(fn):
    stored = mp.mpf(json.load(open(fn))["value"])
    print(f"W_4({k}): mine(fast pieces) - author's dps120 value = {mp.nstr(Wfast.get(3, 0) - stored, 5)}   (mp.dps {mp.mp.dps})")

if slow_dir and os.path.exists(os.path.join(slow_dir, f"slow_red_{k}.pkl")):
    ser = pickle.load(open(os.path.join(slow_dir, f"slow_red_{k}.pkl"), "rb"))
    CAND = {"0": 0, "a": 1 + sp.sqrt(2), "ma": -1 - sp.sqrt(2), "b": sp.sqrt(2) - 1, "mb": 1 - sp.sqrt(2)}
    CANDN = {nm: sp.N(val, 40) for nm, val in CAND.items()}
    slow = {}
    for key, s in ser.items():
        expr = eval(s, dict(sp.__dict__))
        vs = [x for x in expr.free_symbols]
        assert len(vs) == 1, (key, vs)
        v = vs[0]
        P = dict(real={}, ipos=[], ineg=[], const=mp.mpc(0))
        for term in sp.Add.make_args(expr):
            num, den = sp.fraction(sp.together(term))
            assert not num.has(v), ("numerator depends on v", key, term)
            dp = sp.Poly(den, v)
            assert dp.degree() == 1, ("non-simple term", key, term)
            al, be = dp.all_coeffs()
            r = sp.N(-be / al, 40)
            c = num / al
            name = [nm for nm, rv in CANDN.items() if abs(complex(r - rv)) < 1e-25]
            assert len(name) == 1, ("pole not in {0,+-a,+-b}", key, r)
            cn = sp.N(c, int(dps * 1.2) + 20)
            cv = mp.mpc(mp.mpf(str(sp.re(cn))), mp.mpf(str(sp.im(cn))))
            P["real"][name[0]] = P["real"].get(name[0], 0) + cv
        slow[key] = P
    print("slow pieces: keys", sorted(slow))
    if set(slow) != set(pcs): print("KEY SET DIFFERS: fast", sorted(pcs), "slow", sorted(slow))
    worst = mp.mpf(0)
    for key in sorted(set(slow) | set(pcs)):
        fp = pcs.get(key, dict(real={})); sp_ = slow.get(key, dict(real={}))
        poles = set(fp["real"]) | set(sp_["real"])
        dmax = max([abs(fp["real"].get(p, 0) - sp_["real"].get(p, 0)) for p in poles] + [mp.mpf(0)])
        print(f"   key {key}: poles fast {sorted(fp['real'])} slow {sorted(sp_['real'])}  max|coef diff| = {mp.nstr(dmax, 4)}")
        worst = max(worst, dmax)
    print("worst coefficient difference slow vs fast (all keys):", mp.nstr(worst, 4), " mp.dps", mp.mp.dps)
    degs_s = sorted({i + j for (i, j) in slow})
    fss = [piece_integrand(slow, [key for key in slow if sum(key) == d]) for d in degs_s]
    sl = integrate(fss, X, n)
    for d, val in zip(degs_s, sl):
        print(f"  slow  W_{d+1}({k}) = {mp.nstr(val, 45)}   slow - fast = {mp.nstr(val - Wfast.get(d, 0), 5)}  (mp.dps {mp.mp.dps})")
    print(f"RESIDUAL  sum of slow-reducer weights - direct M_{k} = {mp.nstr(sum(sl) - direct, 5)}  (mp.dps {mp.mp.dps})")
print(f"total time {time.time()-t0:.0f}s")
