
import subprocess
r = subprocess.run(['grep','-n','float(x)','step_03.py'], capture_output=True, text=True)
print(r.stdout)
r = subprocess.run(['sed','-n','24,28p','step_03.py'], capture_output=True, text=True)
print(r.stdout)
