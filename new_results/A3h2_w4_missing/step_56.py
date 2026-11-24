import os
wd = '<project>/research/08_theory/A3h2_w4_missing.run'
for f in sorted(os.listdir(wd)):
    p = os.path.join(wd,f)
    print(f, os.path.getsize(p))
