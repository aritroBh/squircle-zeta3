import os
p3 = "../../07_push/P3_fn_graph/claude/result.md"
print(os.path.exists(p3), os.path.getsize(p3) if os.path.exists(p3) else "")
t3 = open(p3).read()
print(t3[:9000])
