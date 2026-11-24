import os
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
steps = sorted([f for f in os.listdir('.') if f.startswith('step_')])
print("step files:", steps)
print("\nprogress.log:")
print(open('progress.log').read())
