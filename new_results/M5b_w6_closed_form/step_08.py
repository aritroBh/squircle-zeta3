import os
base='<project>/research/08_theory/M5a_next_coefficient'
for f in ['step_hp5.py','step_quad5.py','mk_integrands.py','omega5.json']:
    print('='*30,f)
    print(open(os.path.join(base,f)).read())
