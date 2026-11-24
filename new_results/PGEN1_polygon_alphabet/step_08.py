import os
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
idx = tex.find("subsection{Polygons at the Diamond End}")
# next subsection after idx
import re
ms = [m.start() for m in re.finditer(r"\\subsection\{", tex)]
nxt = [m for m in ms if m > idx][0]
sub = tex[idx:nxt]
print(len(sub))
open("/tmp/diamond_end.tex","w").write(sub)
# print in chunks
print(sub[:9000])
