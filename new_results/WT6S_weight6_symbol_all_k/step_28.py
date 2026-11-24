
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
M5A = os.path.join(T, 'M5a_next_coefficient')
print("===== step_quad5.py =====")
print(open(os.path.join(M5A, 'step_quad5.py')).read()[:3000])
