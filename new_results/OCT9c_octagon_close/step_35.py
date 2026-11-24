import subprocess
out = subprocess.run(['sed','-n','1470,1505p','<project>/paper/ieee/main.tex'],capture_output=True,text=True)
print(out.stdout)
