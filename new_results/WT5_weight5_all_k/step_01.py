import os, subprocess
base='<project>'
os.makedirs(base+'/research/08_theory/WT5_weight5_all_k', exist_ok=True)
os.chdir(base+'/research/08_theory/WT5_weight5_all_k')
print(os.getcwd())
# list the theory dir
for d in ['research/08_theory']:
    print(sorted(os.listdir(os.path.join(base,d))))
