import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
A3H5 = os.path.join(ROOT,'research/08_theory/A3h5_exact_certificate')
for s in ['step_29.py','step_30.py','step_31.py','step_32.py','step_33.py','step_34.py']:
    print('='*40, s)
    print(open(os.path.join(A3H5,s)).read())
