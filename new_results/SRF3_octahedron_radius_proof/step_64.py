
import os
base = "<project>"
p = os.path.join(base,"research/08_theory/SRF3R_referee_octahedron/result.md")
txt = open(p).read()
# print the section between the zeros MINOR FIX and "## 4"
i = txt.find("§2/§4.1 zeros")
j = txt.find("## 4")
print(txt[i:j if j>0 else i+9000])
