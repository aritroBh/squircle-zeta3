t = r'''
import re
txt = open("<project>/paper/ieee/main.tex").read()
# find definition of Phi and Psi
for m in re.finditer(r"\\Phi\(c\)", txt):
    s = txt[m.start()-1200:m.start()+300]
    if "Phi" in s and ("=" in s):
        print(s); print("=================")
        break
i = txt.find("Phi(c)")
print(txt[i-2000:i+200])
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout[:7000])
