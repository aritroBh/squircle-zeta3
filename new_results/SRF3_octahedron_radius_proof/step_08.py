
import os, re
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
# find sec:diamond
i = tex.find("sec:diamond}")
print(tex[i-200:i+9000])
