t = r'''
txt = open("<project>/paper/ieee/main.tex").read()
i = txt.find("\\label{tab:closed}")
print(txt[i:i+7000])
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout)
