
import os
base = "<project>"
tow4 = base + "/research/08_theory/TOW4_tower_depth8"
for f in ["tow2lib.py", "tow3lib.py"]:
    print("="*30, f, "="*30)
    print(open(tow4 + "/" + f).read())
