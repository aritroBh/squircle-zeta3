
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
W5 = os.path.join(T, 'WT5_weight5_all_k')
# find wt5_E_pairings
for f in sorted(os.listdir(W5)):
    if 'E_pairings' in f or 'A_kernels' in f or 'A2' in f:
        print(f)
print("===== wt5_E_pairings.py =====")
print(open(os.path.join(W5, 'wt5_E_pairings.py')).read())
