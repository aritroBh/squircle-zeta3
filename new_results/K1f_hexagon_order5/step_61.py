import os
wd = '<project>/research/08_theory/K1f_hexagon_order5'
if os.path.exists(f'{wd}/step_26.out'):
    print(open(f'{wd}/step_26.out').read())
else:
    print("still running")
