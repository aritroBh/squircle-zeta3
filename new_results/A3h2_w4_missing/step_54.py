import os, subprocess
wd = '<project>/research/08_theory/A3h2_missing29'
print(os.path.exists(wd))
if os.path.exists(wd):
    for f in sorted(os.listdir(wd)):
        p = os.path.join(wd,f)
        print(f, os.path.getsize(p))
