
import os
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
print(tex[153688-1500:157500])
