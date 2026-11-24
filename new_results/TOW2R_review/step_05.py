import os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
for fn in ['tow2lib.py']:
    with open(os.path.join(tow2, fn)) as f:
        print('='*30, fn, '='*30)
        print(f.read())
