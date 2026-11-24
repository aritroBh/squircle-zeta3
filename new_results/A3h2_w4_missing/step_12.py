base='<project>'
a3f4=base+'/research/08_theory/A3f4_glanois_basis'
import subprocess
print(subprocess.run(['grep','-l','vals_hi',a3f4+'/step_21.py',a3f4+'/step_22.py'],capture_output=True,text=True).stdout)
# find the file that mentions prec3600
print(subprocess.run(['grep','-rl','3600',a3f4],capture_output=True,text=True).stdout)
