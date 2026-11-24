import subprocess
out = subprocess.run(['grep','-n','prop:order5','<project>/paper/ieee/main.tex'],capture_output=True,text=True)
print(out.stdout)
