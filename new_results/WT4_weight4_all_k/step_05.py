import os
base='<project>/research/08_theory/WT4_weight4_all_k'
for f in ['brief.md','step_01.py','step_02.py','step_03.py']:
    print('='*20, f)
    print(open(os.path.join(base,f)).read()[:5000])
