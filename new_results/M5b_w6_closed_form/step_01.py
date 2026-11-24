import os, subprocess
base='<project>'
print(os.listdir(base))
print(subprocess.run(['ls','-la',os.path.join(base,'research/08_theory/M5a_next_coefficient')],capture_output=True,text=True).stdout)
