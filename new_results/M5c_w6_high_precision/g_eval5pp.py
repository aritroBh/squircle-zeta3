"""Numeric evaluation of gpl_w5pp.tsv vs quadrature (adapted from M5a g_eval5.py)."""
import sys, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
from g_eval_lib import evaluate, mp
t0 = time.time()
val, n_reg = evaluate("gpl_w5pp.tsv")
q = json.load(open("quad5_100.json"))
expected = mp.mpf(q["W5pp"])
got = mp.mpf(val.real.mid().str(110, radius=False, more=True))
err = abs(got - expected)
print(f"gpl_w5pp.tsv: ({time.time()-t0:.0f}s, {n_reg} endpoint-regularized words)")
print(f"   G-value = {mp.nstr(got, 60)}")
print(f"   quad    = {mp.nstr(expected, 60)}")
print(f"   |G - quad| = {mp.nstr(err, 4)}   flint rad = {val.real.rad().str(5)}", flush=True)
json.dump({"G_value": mp.nstr(got, 110), "err_vs_quad": mp.nstr(err, 6),
           "flint_rad": val.real.rad().str(8)}, open("gpl_w5pp_check.json", "w"))
