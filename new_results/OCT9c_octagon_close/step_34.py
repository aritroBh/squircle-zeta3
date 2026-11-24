import subprocess
out = subprocess.run(['grep','-n','edgelaw','<project>/paper/ieee/main.tex'],capture_output=True,text=True)
print(out.stdout[:3000])
