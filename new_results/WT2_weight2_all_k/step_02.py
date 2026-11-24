import os
base = '<project>'
# Create a work directory
work = os.path.join(base, 'wt2_work')
os.makedirs(work, exist_ok=True)
os.chdir(work)
print("workdir:", os.getcwd())

# Read paper lines 836-870
with open(os.path.join(base, 'paper/ieee/main.tex')) as f:
    lines = f.readlines()
print(''.join(lines[830:875]))
