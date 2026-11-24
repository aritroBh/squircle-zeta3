import os
wd = '<project>/research/08_theory/K1f_hexagon_order5'
for f in sorted(os.listdir(wd)):
    p = os.path.join(wd, f)
    print(f, os.path.getsize(p))
