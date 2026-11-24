import os
p = '<project>/research/08_theory/M4c_w5_closed_form/planner'
print(sorted(os.listdir(p)))
print("="*70)
print(open(os.path.join(p,'planner_assemble_M4.py')).read())
