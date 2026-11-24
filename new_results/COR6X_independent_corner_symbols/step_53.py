t = r'''
import re
txt = open("<project>/paper/ieee/main.tex").read()
i = txt.find("sec:corners")
# find the closed-form table
j = txt.find("tab:closed")
print(txt[j-200:j+200])
# search for Li_4 combinations near 2\pi/3
for m in re.finditer(r"Li_4", txt):
    pass
# find eq:R and friends in the closed table
k = txt.find("Table~\\ref{tab:closed}")
print(txt[k-500:k+500] if k>0 else "no ref")
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout)
