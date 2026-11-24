import os
tow4 = '<project>/research/08_theory/TOW4_tower_depth8'
for f in sorted(os.listdir(tow4)):
    p = os.path.join(tow4, f)
    print(f, os.path.getsize(p))
