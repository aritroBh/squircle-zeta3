import os
t = '<project>/research/08_theory/D13t_review_fixes'
f = '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed'
for fn in [ (t,'step_38.py'), (t,'step_37.py'), (t,'step_39.py'), (f,'verify_mS.py') ]:
    print("="*100); print(fn[1], " from ", os.path.basename(fn[0])); print("="*100)
    print(open(os.path.join(*fn)).read())
