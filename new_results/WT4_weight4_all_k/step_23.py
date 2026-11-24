import os
base='<project>/research/08_theory/WTS_weight_k_structure'
for root,dirs,files in os.walk(base):
    depth=root[len(base):].count(os.sep)
    if depth<2:
        print(root, files[:12])
