import os
base = '<project>'
for name in ['ibp_reduce.py']:
    with open(os.path.join(base, 'paper/squircle_note/scripts', name)) as f:
        print("="*30, name, "="*30)
        print(f.read())
