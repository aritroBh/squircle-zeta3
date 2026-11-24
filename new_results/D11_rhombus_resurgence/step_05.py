import subprocess
print(subprocess.run(['python','rhombus_formula.py'],capture_output=True,text=True).stdout)
print(subprocess.run(['python','rhombus_check2.py'],capture_output=True,text=True).stdout)
