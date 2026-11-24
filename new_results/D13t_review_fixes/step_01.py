import os, subprocess
base='research/08_theory/'
for root,dirs,files in os.walk(base):
    for f in files:
        p=os.path.join(root,f)
        print(p, os.path.getsize(p))
