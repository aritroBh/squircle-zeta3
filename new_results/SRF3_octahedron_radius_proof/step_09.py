
import os
base = "<project>"
tex = open(os.path.join(base, "paper/ieee/main.tex")).read()
i = tex.find("The Diamond End")
# find first occurrence of "logit"
j = tex.find("logit")
print("logit idx:", j)
for k in [m.start() for m in __import__('re').finditer(r"logit", tex)]:
    print(k)
