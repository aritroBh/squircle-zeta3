import os
base='research/08_theory/D13s_fix_recertify/'
for root,dirs,files in os.walk(base):
    for f in sorted(files):
        p=os.path.join(root,f)
        print(p, os.path.getsize(p))
