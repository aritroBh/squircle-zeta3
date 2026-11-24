"""verify_const.py <dps> <which> <half>     which in {c0,c1,c2,delta,k3,I:<key>:<pole>}
Computes, with MY OWN integrand evaluation (variable t resp. s=1-t, stable factorised differences, verify_lib) and mpmath tanh-sinh on
[0,1/2], the integral over ONE HALF of the weight-4 integrand combination:
   c_j   : coefficient of k^j in W_4(k) (k>=4), built from my own exact table (verify_lib.exact_table)
   k3    : W_4(3) pieces (red_3.pkl)
   delta : (c0 + 3 c1 + 9 c2) - W_4(3)
half 0: v in [b,1] (variable t),  half 1: v in [1,a] (variable s = 1-t).   Saves $VERIFY_OUT/const_<which>_h<half>_dps<dps>.json (resumable).
PRECISION: mp.mp.dps set after all imports and printed."""
import sys, os, json, time
sys.dont_write_bytecode = True
from verify_lib import *

dps, which, half = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
out_dir = os.environ["VERIFY_OUT"]
fn = os.path.join(out_dir, f"const_{which.replace(':','_')}_h{half}_dps{dps}.json")
if os.path.exists(fn):
    print("exists", fn); sys.exit(0)
set_dps(dps)
print("mp.mp.dps =", mp.mp.dps, "which", which, "half", half, flush=True)

TAB = exact_table()
K3 = k3_table()


def coefs_for(which):
    if which in ("c0", "c1", "c2"):
        j = int(which[1])
        return {kp: v[j] for kp, v in TAB.items()}
    if which == "k3":
        return dict(K3)
    if which == "delta":
        out = {}
        for kp, v in TAB.items():
            p0 = v[0][0] + 3 * v[1][0] + 9 * v[2][0]
            p1 = v[0][1] + 3 * v[1][1] + 9 * v[2][1]
            out[kp] = (p0, p1)
        for kp, (q0, q1) in K3.items():
            a0, a1 = out.get(kp, (Fr(0), Fr(0)))
            out[kp] = (a0 - q0, a1 - q1)
        return out
    if which.startswith("I:"):
        _, ks, pole = which.split(":")
        key = (3, 0) if ks == "K" else (2, 1)
        return {(key, pole): (Fr(1), Fr(0))}
    raise ValueError(which)


pcs = piece_dict(coefs_for(which))
keys = sorted(pcs)
r2, a, b, lam = consts()


def f(t):
    if half == 0:
        v = (2 * t - 1) + mp.sqrt(2 - 4 * t + 4 * t * t)
        l1 = mp.log(t); l2 = mp.log1p(-t)
        dif = {"0": v, "a": v - a, "ma": v + a, "b": 4 * v * t / (v + a), "mb": v + b}
    else:
        w = (2 * t - 1) + mp.sqrt(2 - 4 * t + 4 * t * t)
        v = 1 / w
        l2 = mp.log(t); l1 = mp.log1p(-t)
        dif = {"0": v, "a": -4 * v * t / (v + b), "ma": v + a, "b": v - b, "mb": v + b}
    nd = dict(v=v, dif=dif, l1=l1, l2=l2)
    tot = 0
    for (i, j) in keys:
        tot += eval_pf(pcs[(i, j)], nd) * l1 ** i * l2 ** j
    return tot * 4 * v * v / (v * v + 1)           # dv/dt


t0 = time.time()
val, err = mp.quad(f, [0, mp.mpf(1) / 2], error=True)
print(f"half {half} value = {mp.nstr(val, 40)}   quad err {mp.nstr(err, 5)}   [{time.time()-t0:.0f}s] mp.dps {mp.mp.dps}", flush=True)
json.dump(dict(dps=dps, mp_dps=mp.mp.dps, which=which, half=half, value=mp.nstr(val, dps + 10), err=mp.nstr(err, 5),
               seconds=time.time() - t0), open(fn, "w"))
