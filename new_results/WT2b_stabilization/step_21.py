import os
base='<project>/research/08_theory/WT2_weight2_all_k'
for f in sorted(os.listdir(base)):
    p=os.path.join(base,f)
    print(f, os.path.getsize(p), 'dir' if os.path.isdir(p) else '')
