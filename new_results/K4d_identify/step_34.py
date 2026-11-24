import os
WORK='research/08_theory/K4d_identify'
for lg in ['bgA.log','bgB.log']:
    p=os.path.join(WORK,lg)
    print("="*50, lg)
    if os.path.exists(p):
        print(open(p).read())
    else:
        print("(not yet)")
