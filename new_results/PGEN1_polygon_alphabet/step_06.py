import os, re
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
print(len(tex))
idx = tex.find("Polygons at the Diamond End")
print(idx)
