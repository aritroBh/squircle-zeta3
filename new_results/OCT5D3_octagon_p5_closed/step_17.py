import os
base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
for f in sorted(os.listdir(base)):
    if f.startswith("step") and f.endswith(".py"):
        src = open(base+f).read()
        if 'li32' in src or 'symbol_fit2' in src or 'Li32' in src:
            print("#### ", f)
