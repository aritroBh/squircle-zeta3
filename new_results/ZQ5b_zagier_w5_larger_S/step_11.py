import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
print('=== d21 extra atoms out ===')
print(open('../ZQ5_zagier_weight5/planner/planner_d21_extra_atoms_rels.out').read()[:2000])
print('=== rerun_lll_d17.out ===')
print(open('../ZQ5_zagier_weight5/planner/rerun_lll_d17.out').read()[:1500])
print('=== d21_remaining ===')
print(open('../ZQ5_zagier_weight5/planner/d21_remaining.txt').read()[:1500])
