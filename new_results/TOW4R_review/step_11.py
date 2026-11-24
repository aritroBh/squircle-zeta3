
import os
base = "<project>"
tow2 = base + "/research/08_theory/TOW2_tower_all_depths"
files = sorted(os.listdir(tow2))
for x in files:
    print(f"{os.path.getsize(os.path.join(tow2,x)):>9}  {x}")
