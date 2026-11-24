import os, shutil, time
base = '<project>/research/08_theory'
r = os.path.join(base,'ZQ5_zagier_weight5','out')
print("ZQ5 out:", os.listdir(r) if os.path.exists(r) else "MISSING")
for dd in os.listdir(r):
    p = os.path.join(r,dd)
    print(dd, os.listdir(p) if os.path.isdir(p) else '')
