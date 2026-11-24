import os
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
for fn in ['make_tables.py','gen_functions.py']:
    with open(os.path.join(theory, fn)) as f:
        print('='*30, fn)
        print(f.read()[:5500])
