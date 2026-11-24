import os, subprocess
base = '<project>/research/08_theory'
out = subprocess.run(['find', base + '/CX8_tow2_review', '-name', '*.py'], capture_output=True, text=True)
print(out.stdout)
# also how TOW3 checked the +-i poles for weight 5
tow3 = base + '/TOW3_tower_depth6'
out2 = subprocess.run(['grep','-rln','ipos\|ineg','--include=step_*.py',tow3], capture_output=True, text=True)
print("TOW3 ipos steps:", out2.stdout)
