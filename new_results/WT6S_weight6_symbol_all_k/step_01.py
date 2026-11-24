
import os, subprocess
# First, create the work directory
workdir = 'research/08_theory/WT6S_weight6_symbol_all_k'
os.makedirs(workdir, exist_ok=True)
print("Workdir:", os.path.abspath(workdir))
# List the theory directory
for d in sorted(os.listdir('research/08_theory')):
    print(d)
