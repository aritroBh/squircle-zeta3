import re
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
print(len(tex))
# Find the subsection
idx = tex.find("Polygons at the Diamond End")
print(idx)
