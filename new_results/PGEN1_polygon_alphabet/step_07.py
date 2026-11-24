import os
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
idx = tex.find("Polygons at the Diamond End")
# find end of subsection: next \subsection or \section after idx
import re
m = re.search(r"\\(sub)?section\{", tex[idx+100:])
end = idx+100+m.start() if m else len(tex)
print(tex[idx-200:end])
