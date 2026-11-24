import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
print(os.listdir('../M4e_w5_lchi'))
print('---')
with open('../ZQ5_zagier_weight5/planner/planner_d21_extra_atoms.py') as f:
    print(f.read()[:6000])
