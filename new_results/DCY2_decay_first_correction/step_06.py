base = "<project>/"
txt = open(base+"paper/ieee/main.tex").read()
import re
i = txt.find("How Fast Do the")
print(txt[i:i+200])
# find varpi
for m in re.finditer(r"varpi", txt):
    pass
idxs = [m.start() for m in re.finditer(r"varpi_1", txt)]
print(len(idxs), idxs[:20])
