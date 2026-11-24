
import os
wd = 'research/08_theory/K2b_corner_phase/'
print(sorted(os.listdir(wd)))
print('===== step_35 head =====')
p = wd+'step_35.py'
if os.path.exists(p):
    print(open(p).read()[:2500])
