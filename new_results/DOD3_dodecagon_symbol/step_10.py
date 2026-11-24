import os, subprocess
base = os.path.dirname(os.getcwd())
d1 = os.path.join(base,'DOD1_dodecagon')
# grep for gamma/delta definitions
out = subprocess.run(['grep','-rn','gamdef\|gamma\['+''.join([])+'\|gam\[', d1, '--include=step_*.py','-l'],capture_output=True,text=True)
print(out.stdout)
out = subprocess.run(['grep','-rln','gam', d1, '--include=step_*.py'],capture_output=True,text=True)
print(out.stdout)
