import os
base = '<project>/research/08_theory'
for d in ['ZQ5c_zagier_w5_capped','ZQ5d_zagier_d13_clean','ZQ5_zagier_weight5']:
    p = os.path.join(base,d)
    print('=====', d)
    for root, dirs, files in os.walk(p):
        rel = os.path.relpath(root, p)
        depth = rel.count(os.sep)
        if depth > 2: 
            dirs[:] = []
            continue
        print(rel, ':', sorted(files)[:25])
