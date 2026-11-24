import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
with open('../ZQ5_zagier_weight5/w5lib.py') as f:
    content = f.read()
print(content[14000:])
print('===PLANNER_LI5===')
with open('../ZQ5_zagier_weight5/planner/planner_li5.py') as f:
    print(f.read())
