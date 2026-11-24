# step_16.py - spot-check 3 cells of the planner's sweep file against certified_y38 reruns
import math, json
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
p = '<scratch>'
data = json.load(open(p))
for idx in [0, 20, 41]:
    n1, n2, Cfile = data[idx]
    res = certified_y38(ys(n2), ys(n1), J=24, nsub=4)
    Cre = float(res['C'].upper())
    print(f"cell [{n1:.4g},{n2:.4g}]: file C={Cfile:.4f}, rerun C={Cre:.4f}, match={abs(Cre-Cfile)<1e-6*Cfile}")
with open('progress.log','a') as f:
    f.write("turn4: spot-checked sweep cells idx 0,20,41 - see stdout.\n")
