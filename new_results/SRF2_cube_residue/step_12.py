base = "<project>"
with open(base + "/paper/ieee/main.tex") as f:
    tex = f.read()
print(tex[134998:146667])
