import os
base = "<project>"
with open(os.path.join(base, "paper/ieee/main.tex")) as f:
    tex = f.read()
for needle in ["Ahmed", "Open Problem", "Graphs of"]:
    i = 0
    idxs = []
    while True:
        j = tex.find(needle, i)
        if j < 0: break
        idxs.append(j)
        i = j+1
    print(needle, "->", idxs[:20])
