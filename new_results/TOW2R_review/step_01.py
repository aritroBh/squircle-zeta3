import os, subprocess

# Set up work directory
base = '<project>/research/08_theory'
workdir = os.path.join(base, 'TOW2R_review')
os.makedirs(workdir, exist_ok=True)
print("Workdir:", workdir)

# List the TOW2 directory
tow2 = os.path.join(base, 'TOW2_tower_all_depths')
for root, dirs, files in os.walk(tow2):
    level = root.replace(tow2, '').count(os.sep)
    indent = ' ' * 2 * level
    print(f'{indent}{os.path.basename(root)}/')
    subindent = ' ' * 2 * (level + 1)
    for f in files:
        fp = os.path.join(root, f)
        print(f'{subindent}{f} ({os.path.getsize(fp)} bytes)')
