import subprocess
# Find the paper paragraph
out = subprocess.run(['grep','-n','sign pattern','../../../paper/ieee/main.tex'],capture_output=True,text=True)
print(out.stdout, out.stderr)
out = subprocess.run(['grep','-n','thm:zero','../../../paper/ieee/main.tex'],capture_output=True,text=True)
print(out.stdout[:2000])
