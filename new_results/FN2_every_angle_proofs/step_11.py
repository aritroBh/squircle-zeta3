import os
base = "<project>"
with open(os.path.join(base, "research/07_push/P3_fn_graph/claude/result.md")) as f:
    p3 = f.read()
i = p3.find("Theorem 4")
i0 = p3.rfind("---", 0, i)
i1 = p3.find("## 7", i)
print(p3[i0:i1])
