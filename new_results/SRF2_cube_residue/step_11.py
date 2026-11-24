base = "<project>"
with open(base + "/paper/ieee/main.tex") as f:
    tex = f.read()
# find the section with superellipsoid (3D)
import re
for m in re.finditer(r"\\section\{[^}]*\}", tex):
    print(m.start(), m.group(0))
