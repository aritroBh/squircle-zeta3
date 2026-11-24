import os
base = "<project>"
with open(os.path.join(base, "research/07_push/P3_fn_graph/claude/result.md")) as f:
    p3 = f.read()
print("P3 result.md length:", len(p3))
print(p3[:10000])
