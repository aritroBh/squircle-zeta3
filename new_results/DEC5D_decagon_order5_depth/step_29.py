
import json, os
base = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
# find which step made itab110
import subprocess
for f in sorted(os.listdir(base)):
    if f.startswith("step_") and f.endswith(".py"):
        src = open(base+f).read()
        if "itab" in src:
            print("====", f, "====")
            print(src[:2500])
            print()
