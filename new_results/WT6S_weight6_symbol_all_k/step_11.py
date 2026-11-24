
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
W5 = os.path.join(T, 'WT5_weight5_all_k')
print(open(os.path.join(W5, 'step_50.py')).read())
