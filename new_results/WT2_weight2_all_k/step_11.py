import os
base = '<project>'
for name in ['k5_weight2_rerun.py']:
    with open(os.path.join(base, 'paper/ieee', name)) as f:
        print("="*30, name, "="*30)
        print(f.read())
