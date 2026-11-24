import os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
for fn in ['step_47.py', 'step_48.py', 'step_49.py']:
    with open(os.path.join(tow2, fn)) as f:
        print('='*30, fn, '='*30)
        print(f.read())
