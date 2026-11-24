import time, json
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
import reflib13 as R13

# (a) my own enumeration of maps {0,1/3,1,-1,inf} -> mu_6
found = R13.enum_maps13()
log(f"my exact enumeration over Q(zeta_6): {len(found)} maps")
assert len(found) == 12
for t in found: print("   ", t)

# compare with stored mobius_maps_13.json
import os
base = '<project>/research/08_theory'
stored = sorted(tuple(m['exps']) for m in json.load(open(os.path.join(base,'COR13C_third_corner_certificate','mobius_maps_13.json'))))
log(f"stored maps: {len(stored)}; identical to mine: {stored == found}")
assert stored == found

# (b) task map check: map with exps (0,1,2,4,3) must have f(0)=1, f(1/3)=zeta, f(1)=zeta^2, f(-1)=zeta^4, f(inf)=zeta^3=-1
task = (0,1,2,4,3)
assert task in found
i0 = found.index(task)
log(f"task map exps (0,1,2,4,3): 0->zeta^0=1, 1/3->zeta^1, 1->zeta^2, -1->zeta^4, inf->zeta^3=-1  OK")

# verify directly from the map formula: build map from f(0)=1, f(1)=z^2, f(inf)=z^3 and check f(1/3), f(-1)
MU6 = R13.MU6
M = R13.mobius_from_3(MU6[0], MU6[2], MU6[3])
from fractions import Fraction as Fr
assert R13.qmobius(M, Fr(1,3)) == MU6[1]
assert R13.qmobius(M, Fr(-1)) == MU6[4]
log("direct check: map with f(0)=1,f(1)=z^2,f(inf)=-1 sends 1/3->z, -1->z^4  (exact)")

# (c) omitted-root structure
from collections import Counter
omit = Counter(tuple(sorted(set(range(6)) - set(t)))[0] for t in found)
log(f"omitted-root distribution: {dict(sorted(omit.items()))} (each root omitted exactly twice)")
assert all(v == 2 for v in omit.values())

# (d) shift structure: all maps are shifts (multiply all images by zeta^k) of base tuples
base1 = (0, 1, 2, 4, 3); base2 = (0, 5, 4, 2, 3)
def shifts(t):
    return {tuple((x+k) % 6 for x in t) for k in range(6)}
orbit = shifts(base1) | shifts(base2)
log(f"shift orbits of (0,1,2,4,3) and (0,5,4,2,3) cover all 12 maps: {orbit == set(found)}")
assert orbit == set(found)
with open('progress.log','a') as f:
    f.write("referee step(maps N=6): my exact enumeration = 12 maps == mobius_maps_13.json; task map images verified; structure OK\n")
print("STEP PASS")
