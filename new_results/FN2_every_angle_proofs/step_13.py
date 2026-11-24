import os
base = "<project>"
with open(os.path.join(base, "research/07_push/P3_fn_graph/claude/result.md")) as f:
    p3 = f.read()
print(p3[13500:19000])
