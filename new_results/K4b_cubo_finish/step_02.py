import os
d = "<home>)"
for root, dirs, files in os.walk(d):
    # skip .venv and .git for brevity
    if ".venv" in root or ".git" in root:
        continue
    for f in files:
        if f in ("prev_result.md", "k4_engine.py", "k4_quad.py") or "k4" in f.lower():
            print(os.path.join(root, f))
