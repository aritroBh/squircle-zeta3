import os
files = sorted(os.listdir('.'))
steps = [f for f in files if f.startswith('step_')]
print("deliverables:", [f for f in files if f in ('verify_cor3.py', 'cor3_common.py')])
print("step files:", steps)
print("\nverify_cor3.py header:")
print(open('verify_cor3.py').read()[:400])
