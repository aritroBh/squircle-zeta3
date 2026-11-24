import subprocess
print(subprocess.run(['ls','-la','../DCY2_decay_first_correction/'],capture_output=True,text=True).stdout)
print(subprocess.run(['cat','step_01.py'],capture_output=True,text=True).stdout)
