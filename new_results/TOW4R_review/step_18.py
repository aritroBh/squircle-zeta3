
import os
base = "<project>"
tow4 = base + "/research/08_theory/TOW4_tower_depth8"
for f in ["prove_adj78_fast.py", "prove_adj78.py", "prove_adj78_all.py", "prove_E3b78.py"]:
    print("="*30, f, "="*30)
    print(open(tow4 + "/" + f).read())
