import os
base = '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed'
lines = open(os.path.join(base,'step_38.py')).read().split('\n')
for i, l in enumerate(lines):
    if 85 <= i <= 150:
        print(f"{i+1:4d}: {l}")
