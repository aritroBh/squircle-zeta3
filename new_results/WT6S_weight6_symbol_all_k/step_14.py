
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
M5A = os.path.join(T, 'M5a_next_coefficient')
for f in ['w6_atoms.py', 'w6_solve.py']:
    print(f"===== {f} =====")
    print(open(os.path.join(M5A, f)).read())
    print()
