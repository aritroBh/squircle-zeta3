import os, subprocess

# Set up work directory
workdir = '<project>/research/08_theory/TOW5_tower_depth10'
os.makedirs(workdir, exist_ok=True)
os.chdir(workdir)
print("CWD:", os.getcwd())

# Initialize progress.log
with open('progress.log','w') as f:
    f.write("# TOW5 progress log\n")

# List the theory directory to see what's there
base = '<project>/research/08_theory'
for d in sorted(os.listdir(base)):
    print(d)
