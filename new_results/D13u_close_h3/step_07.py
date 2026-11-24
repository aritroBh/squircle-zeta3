import os
base = '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed'
for f in ['step_18.py','certify_H3.py']:
    print("="*100)
    print(f)
    print("="*100)
    print(open(os.path.join(base,f)).read())
