import os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
with open(os.path.join(tow2, 'genintegrand.py')) as f:
    print(f.read())
with open(os.path.join(tow2, 'progress.log')) as f:
    print('=== TOW2 progress.log ===')
    print(f.read())
