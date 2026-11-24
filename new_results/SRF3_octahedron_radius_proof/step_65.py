
import os
base = "<project>"
txt = open(os.path.join(base,"research/08_theory/SRF3R_referee_octahedron/result.md")).read()
i = txt.find("## 4")
j = txt.find("## 5")
print(txt[i:j])
