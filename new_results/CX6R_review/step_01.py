import os, subprocess
print(os.getcwd())
print(subprocess.run(['ls','-la'],capture_output=True,text=True).stdout)
print(subprocess.run(['ls','-la','<project>/research/08_theory/CX6_decay_rest/'],capture_output=True,text=True).stdout)
