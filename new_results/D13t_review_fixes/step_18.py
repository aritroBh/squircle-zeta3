import difflib, subprocess
base='research/08_theory/D13s_fix_recertify/D13s_work/'
for pair in [('buggy/step_38.py','fixed/step_38.py'),
             ('D13h_final_cell__step_10_J192.py','fixed/step_10_J192.py'),
             ('buggy/step_18.py','fixed/step_18.py')]:
    a=open(base+pair[0]).read().split('\n'); b=open(base+pair[1]).read().split('\n')
    print('='*25, pair)
    for l in difflib.unified_diff(a,b,lineterm='',n=1):
        if l.startswith(('+','-')) and not l.startswith(('+++','---')):
            print(l[:200])
