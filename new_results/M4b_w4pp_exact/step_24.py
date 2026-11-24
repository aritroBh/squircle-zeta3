import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
A3H5 = os.path.join(ROOT,'research/08_theory/A3h5_exact_certificate')
print(open(os.path.join(A3H5,'step_32.py')).read())
