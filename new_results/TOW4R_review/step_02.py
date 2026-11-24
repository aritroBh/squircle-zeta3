
import os
base = "<project>"
tow4 = os.path.join(base, "research/08_theory/TOW4_tower_depth8")
for x in sorted(os.listdir(tow4)):
    full = os.path.join(tow4, x)
    print(f"{os.path.getsize(full):>10}  {x}")
