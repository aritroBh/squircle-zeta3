import os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
for fn in ['step_51.py', 'step_52.py', 'step_55.py', 'step_58.py', 'step_43.py']:
    with open(os.path.join(tow2, fn)) as f:
        print('='*30, fn, '='*30)
        print(f.read())
