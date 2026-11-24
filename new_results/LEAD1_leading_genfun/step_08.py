import os
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
for fn in ['check_top.py']:
    with open(os.path.join(theory, fn)) as f:
        print('='*30, fn)
        print(f.read())
