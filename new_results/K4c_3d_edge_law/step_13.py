
import subprocess
r = subprocess.run(['grep','-n','np.array(v','step_03.py'], capture_output=True, text=True)
print(r.stdout)
r = subprocess.run(['sed','-n','18,22p','step_03.py'], capture_output=True, text=True)
print(r.stdout)
