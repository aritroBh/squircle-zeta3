
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
W5 = os.path.join(T, 'WT5_weight5_all_k')
for f in ['w5_cobracket.py', 'cobracket_e.json']:
    print(f"===== {f} =====")
    print(open(os.path.join(W5, f)).read())
    print()
