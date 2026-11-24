import os
wts = '<project>/research/08_theory/WTS_weight_k_structure'
for root, dirs, files in os.walk(wts):
    depth = root.replace(wts,'').count('/')
    if depth > 1: continue
    for f in sorted(files):
        p = os.path.join(root,f)
        print(os.path.relpath(p, wts), os.path.getsize(p))
