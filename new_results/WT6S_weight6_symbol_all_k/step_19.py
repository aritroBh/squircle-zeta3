
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
W5 = os.path.join(T, 'WT5_weight5_all_k')
print("===== deltas_deg4.json =====")
print(open(os.path.join(W5, 'deltas_deg4.json')).read()[:2500])
print("\n===== wt5_A_kernels.py =====")
print(open(os.path.join(W5, 'wt5_A_kernels.py')).read()[:4000])
