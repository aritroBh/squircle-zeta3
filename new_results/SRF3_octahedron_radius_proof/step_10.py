
import os
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
print(tex[62000:75000])
