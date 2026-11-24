
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
M5A = os.path.join(T, 'M5a_next_coefficient')
print("===== gexpr5.py =====")
print(open(os.path.join(M5A, 'gexpr5.py')).read()[:2500])
