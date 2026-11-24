import os
base = "<project>"
with open(os.path.join(base, "research/07_push/P3_fn_graph/claude/result.md")) as f:
    p3 = f.read()
# Print the truncated middle part: Prop 2 statement, and Theorem 4
i2 = p3.find("Proposition 2")
print(p3[3000:10000])
