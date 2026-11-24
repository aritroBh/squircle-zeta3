import os, subprocess
base = '<project>/research/08_theory'
# List relevant directories
for d in ['D13t_review_fixes', 'D13s_fix_recertify/D13s_work/fixed']:
    p = os.path.join(base, d)
    print("=== ", p)
    for f in os.listdir(p):
        print("   ", f, os.path.getsize(os.path.join(p,f)))
