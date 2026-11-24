import os
base = '<project>'
for name in ['mk_integrands.py', 'qsqrt2.py']:
    with open(os.path.join(base, 'paper/squircle_note/scripts', name)) as f:
        print("="*30, name, "="*30)
        print(f.read())
