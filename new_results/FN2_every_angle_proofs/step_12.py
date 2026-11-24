import os
base = "<project>"
with open(os.path.join(base, "research/07_push/P3_fn_graph/claude/result.md")) as f:
    p3 = f.read()
for needle in ["Theorem 4", "## 6", "Lambert", "contact", "Contact"]:
    idxs = []
    st = 0
    while True:
        j = p3.find(needle, st)
        if j<0: break
        idxs.append(j); st=j+1
    print(needle, idxs[:15])
