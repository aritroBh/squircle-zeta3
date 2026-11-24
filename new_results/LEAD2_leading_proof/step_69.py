import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
files = sorted(os.listdir('.'))
steps = sorted([f for f in files if f.startswith('step_') and f.endswith('.py')])
print(len(steps), 'step files:', steps)
print()
print(open('progress.log').read())
