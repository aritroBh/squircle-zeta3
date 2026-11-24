import subprocess, os
d = 'research/08_theory/A3f3_level8_compute'
out = subprocess.run(['grep','-rl','mpmath',d],capture_output=True,text=True).stdout
print(out)
out2 = subprocess.run(['grep','-rn','freevals',d],capture_output=True,text=True).stdout
print(out2[:3000])
