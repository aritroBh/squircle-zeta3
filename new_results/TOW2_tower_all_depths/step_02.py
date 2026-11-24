import os
base = '<project>/research/08_theory'
for d in ['TOW1_tower_law', 'WTS_weight_k_structure', 'LEAD1_leading_genfun', 'TOW2_tower_all_depths']:
    p = os.path.join(base, d)
    print('===', d)
    for root, dirs, files in os.walk(p):
        for f in sorted(files):
            fp = os.path.join(root, f)
            print('  ', os.path.relpath(fp, p), os.path.getsize(fp))
