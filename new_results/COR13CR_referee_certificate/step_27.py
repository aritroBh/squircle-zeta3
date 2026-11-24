import os
base = '<project>/research/08_theory'
for s in ['step_06.py','step_07.py']:
    print("="*40, s)
    print(open(os.path.join(base,'HEXC2R_referee_certificate','steps',s)).read())
