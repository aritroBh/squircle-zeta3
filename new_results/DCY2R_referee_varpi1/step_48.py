import subprocess
r = subprocess.run(['bash','-c','sed -n "1240,1295p" ../../../paper/ieee/main.tex'],capture_output=True,text=True)
print(r.stdout)
