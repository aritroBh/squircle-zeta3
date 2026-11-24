import os
target = '<project>/research/08_theory/COR8B_depth_two_two_terms'
for f in ['step_04.py', 'step_09.py']:
    print('='*30, f)
    with open(os.path.join(target, f)) as fh:
        print(fh.read())
