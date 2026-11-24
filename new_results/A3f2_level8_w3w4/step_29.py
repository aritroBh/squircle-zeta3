import subprocess
print(subprocess.run(['sed','-n','120,220p','<project>/paper/squircle_note/scripts/qsqrt2.py'],capture_output=True,text=True).stdout)
