import os
base = "<project>"
with open(os.path.join(base, "paper/ieee/main.tex")) as f:
    tex = f.read()
print("len:", len(tex))
idx = tex.find("sec:fn")
print("sec:fn label at", idx)
# print surrounding section
start = tex.rfind("\\section", 0, idx)
print("section start:", start, tex[start:start+200])
