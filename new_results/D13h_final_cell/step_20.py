# step_17.py - junction consistency + J-refinement stability of the final cell
import math, subprocess, sys, json
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
# (1) certified_y38 on a box just above 1e10: N in [1e10, 1.05e10]
r = certified_y38(ys(1.05e10), ys(1e10), J=24, nsub=2)
print("junction box [1e10,1.05e10]: step_38 C =", float(r['C'].upper()))
print("  (final cell covers [1e10,oo) with C=57.0676 -- must be >= the limit of these)")

# (2) J-refinement of final cell: patch J=96, zstar=3.0 and rerun as a subprocess variant
src = open('step_10.py').read()
v = src.replace('J = 48; zstar = arb("3.0")', 'J = 96; zstar = arb("3.0")')
v = v.replace('open("final_cell.json","w")', 'open("final_cell_J96.json","w")')
open('step_10_J96.py','w').write(v)
r2 = subprocess.run([sys.executable, 'step_10_J96.py'], capture_output=True, text=True, timeout=300)
for line in r2.stdout.splitlines():
    if 'FINAL CELL' in line or 'W2reg' in line or 'W3reg' in line:
        print(line)
print("stderr:", r2.stderr[-500:] if r2.stderr else "(none)")
