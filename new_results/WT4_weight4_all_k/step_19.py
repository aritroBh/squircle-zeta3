import os
base='<project>/research/08_theory/M4b2_w4pp_nice'
print(open(os.path.join(base,'progress.log')).read()[-3000:])
print("===== m4b2_common.py")
print(open(os.path.join(base,'m4b2_common.py')).read()[:3000])
