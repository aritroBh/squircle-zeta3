import os
# My work directory
work = 'research/08_theory/OCT4_octagon_letters'
os.chdir(work)
print("CWD:", os.getcwd())

# Read existing steps
for f in ['step_01.py', 'step_02.py']:
    print('='*50, f)
    with open(f) as fh:
        print(fh.read())
