t = r'''
import re
txt = open("<project>/paper/ieee/main.tex").read()
# Find Table tab:closed
i = txt.find("\\label{tab:closed}")
print(txt[i-5000:i+100])
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout[:7000])
