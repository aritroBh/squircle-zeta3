import os
base = '<project>/research/08_theory'
for s in ['step_05.py','step_09.py']:
    print("="*40, s)
    src = open(os.path.join(base,'HEXC2R_referee_certificate','steps',s)).read()
    print(src[:4800])
