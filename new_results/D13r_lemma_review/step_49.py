import subprocess, json
d='research/08_theory/D13h_final_cell'
cells = json.load(open(d+'/sweep_cells_400_1e10.json'))
script = """
import sys, math, json
sys.path.insert(0,'.')
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
data = %s
for idx in [0, 20, 41, 42]:
    n1, n2, Cfile = data[idx]
    res = certified_y38(ys(n2), ys(n1), J=24, nsub=4)
    Cre = float(res['C'].upper())
    print(f"cell {idx} [{n1:.6g},{n2:.6g}]: file C={Cfile:.4f}, rerun C={Cre:.4f}, match={abs(Cre-Cfile)<1e-9*Cfile}")
""" % (repr(cells),)
r = subprocess.run(['<home>)/.venv/bin/python','-c',script],
                   cwd=d, capture_output=True, text=True)
print(r.stdout, r.stderr[-300:] if r.stderr else "")
