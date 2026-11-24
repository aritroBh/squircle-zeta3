base='<project>'
W=base+'/research/08_theory/WT5_weight5_all_k/'
script = '''"""WT5 step C2: e1/e2/e3 TSVs vs Arb quadrature (800-bit values) at ~210 digits + W5pp control."""
import sys, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
from g_eval_lib import evaluate, mp
t0 = time.time()
out = {}
for name, tsv, quadjson in [("e1","e1.tsv","e1_800.json"), ("e2","e2.tsv","e2_800.json"), ("e3","e3.tsv","e3_800.json"),
                            ("w5pp","gpl_w5pp.tsv","w5pp_800.json")]:
    val, n_reg = evaluate(tsv)
    q = json.load(open(quadjson))
    expected = mp.mpf(q["value"])
    got = mp.mpf(val.real.mid().str(250, radius=False, more=True))
    err = abs(got - expected)
    print(f"{tsv}: ({time.time()-t0:.0f}s, {n_reg} endpoint-regularized words)")
    print(f"   |G - quad| = {mp.nstr(err, 4)}   flint rad = {val.real.rad().str(5)}", flush=True)
    out[name] = {"G_value": mp.nstr(got, 210), "err_vs_quad": mp.nstr(err, 6), "flint_rad": val.real.rad().str(8)}
json.dump(out, open("check_tsv_vs_quad.json", "w"), indent=1)
with open("progress.log","a") as f:
    f.write("step_C2 check(210d): " + "; ".join(f"{k} |G-quad|={v['err_vs_quad']}" for k,v in out.items()) + "\\n")
print("DONE")
'''
open(W+'wt5_C2_check.py','w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_C2_check.py'], cwd=W, capture_output=True, text=True, timeout=3000)
print(r.stdout); print('STDERR:', r.stderr[-800:])
