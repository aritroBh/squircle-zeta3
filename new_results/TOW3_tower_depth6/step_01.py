import os, subprocess
base = '<project>/research/08_theory'
for d in ['TOW3_tower_depth6', 'TOW2_tower_all_depths', 'TOW1_tower_law', 'SUB1_subleading_table', 'TOWN1_tower_large_n']:
    p = os.path.join(base, d)
    print('='*80)
    print(d, os.path.exists(p))
    if os.path.exists(p):
        for f in sorted(os.listdir(p)):
            print('  ', f)
