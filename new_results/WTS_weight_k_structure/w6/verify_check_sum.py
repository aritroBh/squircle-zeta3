"""verify step: (1) direct-from-definition M_k vs  sum_d W_{d+1}(k)  (own Gauss-Legendre quadrature of the
exact finals_k.json), for given k;  (2) own W_6(k) vs the stored 110-dps value.
usage: verify_check_sum.py TARGET_DIGITS K [K ...]
Everything is computed here; nothing is taken from s02_eval.py (except reading the stored value for comparison).
"""
import sys, os, json, math, time
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
import verify_lib as VL
import mpmath as mp

TGT = int(sys.argv[1])
KS = [int(x) for x in sys.argv[2:]]
Yq = math.log(10) * (TGT + 14) + 30           # quadrature cut for the pieces

for k in KS:
    t0 = time.time()
    # ---- sum of all weight pieces from the exact finals (own quadrature)
    VL.set_dps(TGT + 10)
    fin = VL.load_finals(os.path.join(W, f"finals_{k}.json"))
    Wd = {}
    for d in range(1, k + 1):
        val, _ = VL.weight_part(fin, d, TGT + 10)
        Wd[d] = val
    keys0 = [kk for kk in fin if sum(kk) == 0]
    assert not keys0, "unexpected weight-1 keys"
    S = sum(Wd.values(), mp.mpf(0))
    print(f"k={k}: mp.mp.dps={mp.mp.dps}  pieces done [{time.time()-t0:.0f}s]", flush=True)
    for d in range(1, k + 1):
        print(f"   W_{d+1}({k}) = {mp.nstr(Wd[d], 40)}", flush=True)
    # ---- direct M_k at a working precision that absorbs the t->0 cancellation
    Y = math.log(10) * (TGT + 12) + 25 + 3 * k
    wd = TGT + int(Y / 2.302585) + 30
    VL.set_dps(wd)
    Md = VL.direct_M(k, TGT, Y=Y)
    VL.set_dps(TGT + 10)
    print(f"k={k}: direct M_k (def. P(p)/8 expansion), working mp.mp.dps={wd}, Y={Y:.1f} -> {mp.nstr(Md, 50)}", flush=True)
    print(f"k={k}: sum_d W_(d+1) (mp.mp.dps={mp.mp.dps})                           -> {mp.nstr(S, 50)}", flush=True)
    print(f"k={k}: RESIDUAL direct - sum = {mp.nstr(Md - S, 5)}   [{time.time()-t0:.0f}s]", flush=True)
    # ---- compare W_6(k) with stored 110-dps value
    if 5 <= k:
        fn = os.path.join(W, "vals_110", f"k{k}_d5.json")
        if os.path.exists(fn):
            st = json.load(open(fn))
            stv = mp.mpf(st["value"])
            print(f"k={k}: own W_6 = {mp.nstr(Wd[5], 60)}", flush=True)
            print(f"k={k}: stored W_6 (110-dps run) = {mp.nstr(stv, 60)}", flush=True)
            print(f"k={k}: W_6 own - stored = {mp.nstr(Wd[5] - stv, 5)}  (own mp.mp.dps={mp.mp.dps})", flush=True)
