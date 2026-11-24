import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
print(open('progress.log').read())
print('--- files ---')
for f in sorted(os.listdir('.')):
    if not f.startswith('step_'): print(f, os.path.getsize(f))
