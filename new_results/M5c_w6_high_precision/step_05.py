import os, shutil
base='<project>/research/08_theory/M5b_w6_closed_form'
work='<project>/research/08_theory/M5c_w6_high_precision'
print(open(os.path.join(base,'result.md')).read()[6000:9000])
print("===== progress.log tail =====")
log=open(os.path.join(base,'progress.log')).read().strip().splitlines()
for l in log[-40:]: print(l)
