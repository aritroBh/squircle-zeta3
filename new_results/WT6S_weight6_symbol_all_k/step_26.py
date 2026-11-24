
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
M5A = os.path.join(T, 'M5a_next_coefficient')
print(open(os.path.join(M5A, 'mk_integrands.py')).read()[:3000])
