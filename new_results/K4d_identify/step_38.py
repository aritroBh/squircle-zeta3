import os
WORK='research/08_theory/K4d_identify'
for lg in ['bgA.log','bgB.log']:
    p=os.path.join(WORK,lg)
    print("="*40, lg)
    print(open(p).read() if os.path.exists(p) else "(none)")
# check job outputs
for o in ['jobA_out.txt','jobB_out.txt']:
    p=os.path.join(WORK,o)
    if os.path.exists(p):
        t=open(p).read()
        if t.strip(): print("---",o,"---\n",t[-800:])
