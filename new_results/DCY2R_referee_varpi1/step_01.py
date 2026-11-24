import os, subprocess
print(os.getcwd())
print(subprocess.run(['ls','-la'],capture_output=True,text=True).stdout)
print(subprocess.run(['ls','-la','../'],capture_output=True,text=True).stdout)
