"""Compare quadrature results with the known closed forms.  Usage: python w4_controls.py <dps>"""
import sys, glob, json
from w4_common import *
from w4_known import *
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 90
set_dps(dps)
print("mp.mp.dps =", mp.mp.dps)
def val(w, k, dp=dps):
    fn = os.path.join(RES, f"w{w}_k{k}_dps{dp}.json")
    return mp.mpf(json.load(open(fn))["value"]) if os.path.exists(fn) else None
rows = []
def report(name, num, ref, tol=mp.mpf(10)**(-(dps - 12))):
    if num is None: return
    r = num - ref
    ok = abs(r) < tol
    print(f"{name:34s} residual {mp.nstr(r, 5):>14s}   dps {mp.mp.dps}   {'PASS' if ok else 'FAIL'}")
    rows.append({"name": name, "residual": mp.nstr(r, 5), "pass": bool(ok)})
for k in range(2, 11):
    report(f"W_2({k}) - (pi^2+6 log^2 2)/24", val(2, k), W2_exact())
for k in range(2, 11):
    report(f"W_3({k}) - [W3' - (k-3) S]", val(3, k), W3_known(k))
report("W_4(3) - W4 (thm:W4 closed form)", val(4, 3), W4_thm())
report("W_4(4) - W4'' (W4pp_800.txt)", val(4, 4), W4pp_txt())
report("W_4(4) - W4'' (w4pp_closed_form.json)", val(4, 4), W4pp_closed())
print("sanity: closed form W4'' vs txt:", mp.nstr(W4pp_closed() - W4pp_txt(), 5))
print("sanity: W4 thm vs tex 48-digit value 0.795785059832240811320359819629726762546527413871 :",
      mp.nstr(W4_thm() - mp.mpf("0.795785059832240811320359819629726762546527413871"), 5))
json.dump(rows, open(os.path.join(RES, f"controls_dps{dps}.json"), "w"), indent=1)
