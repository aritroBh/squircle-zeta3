import os
M4E = '<project>/research/08_theory/M4e_w5_lchi'
for f in sorted(os.listdir(M4E)):
    if f.startswith('step_') and f.endswith('.py') and f not in ['step_01.py','step_02.py','step_03.py']:
        print("="*70)
        print(f)
        print(open(os.path.join(M4E,f)).read()[:2500])
