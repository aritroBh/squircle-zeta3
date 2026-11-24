import subprocess, sys
r = subprocess.run(['grep','-n','truediv\\|rtruediv','w5lib.py'],capture_output=True,text=True)
print(r.stdout)
r2 = subprocess.run(['sed','-n','30,42p','w5lib.py'],capture_output=True,text=True)
print(r2.stdout)
