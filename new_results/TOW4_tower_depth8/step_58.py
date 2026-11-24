import os, subprocess
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
out = subprocess.run(['grep','-rln','38939\|tower law\|cells, 0 bad','--include=step_*.py',tow3], capture_output=True, text=True)
print(out.stdout)
out2 = subprocess.run(['grep','-rln','38939', tow3], capture_output=True, text=True)
print("md/pkl:", out2.stdout)
