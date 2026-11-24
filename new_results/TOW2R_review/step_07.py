import os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
for fn in ['prove_S12.py', 'prove_S34.py', 'prove_tops.py']:
    with open(os.path.join(tow2, fn)) as f:
        print('='*30, fn, '='*30)
        print(f.read())
