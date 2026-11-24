p = "<project>/paper/ieee/main.tex"
with open(p) as f:
    tex = f.read()
print(tex[119500:122650])
