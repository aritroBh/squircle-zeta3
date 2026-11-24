import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, time
import mpmath as mp
import w5lib as W
J = json.load(open("../M4i_w5_noL/L_relation.json"))
terms = J["terms"]
P = 830
mp.mp.dps = P+30
lam = mp.log(2); L = mp.log(1+mp.sqrt(2)); pi2 = mp.pi**2
def atom_val(key):
    if key == "z5": return mp.zeta(5)
    if key.startswith("Li5("):
        inner = key[4:-1]
        x = eval(inner, {"sqrt2": mp.sqrt(2)})
        return W.li5_series(x, P)
    v = mp.mpf(1)
    for f in key.split("*"):
        v *= {"lam": lam, "L": L, "pi2": pi2}[f]
    return v
rhs = sum(Fraction(val[1])*atom_val(key) for key, val in terms.items())
Tv = W.Lchi_hurwitz(8, P)/mp.sqrt(2)
print("dps now:", mp.mp.dps)
print("rhs - L/sqrt2 (computed at 860d):", mp.nstr(rhs - Tv, 5))
# also check the worst-conditioned single Li5
x = eval("(1/2+1/4*sqrt2)", {"sqrt2": mp.sqrt(2)})
v1 = W.li5_series(x, P); v2 = mp.polylog(5, x)
print("series-polylog for x=(2+sqrt2)/4:", mp.nstr(v1-v2, 3))
'''
open('steps/step_09_m4eval.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_09_m4eval.py'],capture_output=True,text=True,timeout=1200)
print(r.stdout); print("STDERR:", r.stderr[-1000:])
with open('progress.log','a') as f: f.write("step_09_m4eval: M4i RHS residual redone at 860d\n")
