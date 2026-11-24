import subprocess
print(subprocess.run(['sed','-n','1,120p','<project>/paper/squircle_note/scripts/qsqrt2.py'],capture_output=True,text=True).stdout)
