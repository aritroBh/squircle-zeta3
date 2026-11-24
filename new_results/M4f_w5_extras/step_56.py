import os
T='<project>/research/08_theory'
print(os.listdir(f'{T}/M4c_w5_closed_form/planner'))
print('='*80)
print(open(f'{T}/M4c_w5_closed_form/planner/planner_assemble_M4.py').read()[:5000])
