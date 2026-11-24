p = "<project>/paper/ieee/main.tex"
with open(p) as f:
    tex = f.read()
i = tex.find("Polygons at the Diamond End")
print(tex[i:i+7500])
