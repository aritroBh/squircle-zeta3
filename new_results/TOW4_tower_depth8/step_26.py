import os
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
# find the rho check step
import subprocess
out = subprocess.run(['grep','-l','720', tow3+'/step_'+str(i).zfill(2)+'.py'] if False else ['grep','-rl','720','--include=step_*.py',tow3], capture_output=True, text=True)
print(out.stdout)
out = subprocess.run(['grep','-rn','rho_level\|initial\|k0 =\|k = d\b','--include=step_3*.py',tow3], capture_output=True, text=True)
print(out.stdout[:2000])
