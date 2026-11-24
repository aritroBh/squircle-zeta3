base = "<project>"
import re
with open(base + "/paper/ieee/main.tex") as f:
    tex = f.read()
# find "The Superellipsoid" section and lem:holo
i = tex.find("Superellipsoid")
j = tex.find("lem:holo")
print("idx:", i, j)
# print section around Superellipsoid
k = tex.find("\\section", max(0, i-2000))
print(tex[k:k+200])
