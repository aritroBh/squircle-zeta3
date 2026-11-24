base = "<project>"
with open(base + "/paper/ieee/main.tex") as f:
    tex = f.read()
i = tex.find("The area splits the same way")
print(tex[i-200:i+4500])
