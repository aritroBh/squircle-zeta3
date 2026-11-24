t = r'''
import re
txt = open("<project>/paper/ieee/main.tex").read()
for m in re.finditer(r"thm:W4R", txt):
    s = txt[m.start()-2500:m.start()+1500]
    if "Theorem" in s or "theorem" in s:
        print(s)
        print("=========================")
        break
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout[:6000])
