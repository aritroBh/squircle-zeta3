import os, pickle, sys
base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
print([f for f in os.listdir(base) if 'atom' in f or 'kernel' in f or 'block' in f])
for f in ["step_76.py","step_77.py"]:
    print("="*100); print(f)
    print(open(base+f).read()[:3500])
