
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
M5A = os.path.join(T, 'M5a_next_coefficient')
print("===== g_eval5.py =====")
print(open(os.path.join(M5A, 'g_eval5.py')).read()[:4000])
