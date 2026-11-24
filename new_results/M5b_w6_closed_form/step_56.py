import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
m5a='<project>/research/08_theory/M5a_next_coefficient'
src=open(os.path.join(m5a,'g_eval5.py')).read()
# guard the bottom driver part
old = src[src.index('t0 = time.time()'):]
lib = src.replace(old, '')
open(os.path.join(m5b,'g_eval_lib.py'),'w').write(lib)
print('lib written; tail removed was %d chars' % len(old))
ev = '''"""Numeric evaluation of gpl_w5pp.tsv vs quadrature (adapted from M5a g_eval5.py)."""
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
'''
open(os.path.join(m5b,'g_eval5pp.py'),'w').write(ev)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','g_eval5pp.py'],capture_output=True,text=True,cwd=m5b,timeout=1500)
print(r.stdout[-700:]); print(r.stderr[-400:])
