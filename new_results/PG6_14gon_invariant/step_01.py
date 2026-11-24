import subprocess, os
# Set up working directory
workdir = 'research/08_theory/PG6_14gon_invariant/'
os.makedirs(workdir, exist_ok=True)
print("CWD:", os.getcwd())
# List the theory folder to see what's there
for p in sorted(os.listdir('research/08_theory')):
    print(p)
