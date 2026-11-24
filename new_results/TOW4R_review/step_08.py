
import os
base = "<project>"
tow4 = base + "/research/08_theory/TOW4_tower_depth8"
for f in ["tow4lib.py", "forced78.py"]:
    print("="*30, f, "="*30)
    print(open(tow4 + "/" + f).read())
