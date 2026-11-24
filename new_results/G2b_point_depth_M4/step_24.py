
import subprocess, re
tex = open("<project>/paper/ieee/main.tex").read()
i = tex.find("sec:M4")
print(tex[i-200:i+5000])
