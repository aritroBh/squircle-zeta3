import os
p = '<project>/research/08_theory/M4e_w5_lchi'
for f in ['step_01.py','step_02.py','step_03.py']:
    print("="*70)
    print(f)
    print(open(os.path.join(p,f)).read())
