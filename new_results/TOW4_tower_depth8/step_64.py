import os, subprocess
base = '<project>/research/08_theory'
out = subprocess.run(['find', base, '-name', 'kern_*.pkl'], capture_output=True, text=True)
print(out.stdout)
# how did TOW3 do the k=5 kernel check? find the step
tow3 = base + '/TOW3_tower_depth6'
out2 = subprocess.run(['grep','-rln','kern_','--include=step_*.py',tow3], capture_output=True, text=True)
print(out2.stdout)
