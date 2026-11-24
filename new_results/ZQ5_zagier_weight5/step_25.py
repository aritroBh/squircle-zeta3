import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
import json, re
import mpmath as mp
import w5lib as W
from pathlib import Path
DPS = 830
mp.mp.dps = DPS + 40
# --- exec planner helpers (same as planner_check_Lrel does) ---
srcp = Path("../M4e_w5_lchi/planner/planner_eval_terms.py").read_text()
srcp = srcp[:srcp.index('import math\nterms =')]
import sys as _s
_s.argv = ['x', '../M4i_w5_noL/L_relation.json']
g = {'__name__': 'helpers'}
# their helper file references TERMS/DPS inside; we only need atom(), li(), field(), num(), split(), Lchi
exec(compile(srcp, 'helpers', 'exec'), g)
atom, num, split, Lchi = g['atom'], g['num'], g['split'], g['Lchi']
# --- my evaluation ---
lam = mp.log(2); L = mp.log(1+mp.sqrt(2)); pi2 = mp.pi**2
def my_atom(key):
    if key == "z5": return mp.zeta(5)
    if key.startswith("Li5("):
        inner = key[4:-1]
        x = eval(inner, {"sqrt2": mp.sqrt(2)})
        return W.li5_series(x, DPS)
    v = mp.mpf(1)
    for f in key.split("*"):
        v *= {"lam": lam, "L": L, "pi2": pi2}[f]
    return v
terms = json.load(open("../M4i_w5_noL/L_relation.json"))["terms"]
tot_mine = mp.mpf(0); tot_theirs = mp.mpf(0)
for key, (p, q) in terms.items():
    c = num(F(p), F(q))
    vm = my_atom(key)
    vt = mp.mpf(1)
    for f in split(key): vt *= atom(f)
    d = vm - vt
    tot_mine += c*vm; tot_theirs += c*vt
    flag = " <<<" if abs(d) > mp.mpf('1e-860') else ""
    print(f"{key:>26s}  mine-theirs = {mp.nstr(d,3)}{flag}")
print("total mine - theirs:", mp.nstr(tot_mine - tot_theirs, 3))
print("theirs - Lchi:", mp.nstr(tot_theirs - Lchi, 3))
'''
open('steps/step_12_atomdiff.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_12_atomdiff.py'],capture_output=True,text=True,timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_12_atomdiff: per-atom diff vs planner\n")
