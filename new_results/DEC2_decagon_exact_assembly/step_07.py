import re
p = "<project>/paper/ieee/main.tex"
with open(p) as f:
    tex = f.read()
i = tex.find("Diamond End")
print(tex[i-200:i+5200])
