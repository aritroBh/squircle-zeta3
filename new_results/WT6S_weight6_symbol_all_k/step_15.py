
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
M5A = os.path.join(T, 'M5a_next_coefficient')
print("===== planner6_pairing.py =====")
print(open(os.path.join(M5A, 'planner6_pairing.py')).read())
print("\n===== M5b n6_solution.json =====")
M5B = os.path.join(T, 'M5b_w6_closed_form')
print(open(os.path.join(M5B, 'n6_solution.json')).read())
