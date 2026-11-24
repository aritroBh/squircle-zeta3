import os
d = "<home>)"
print(os.listdir(d))
for f in ["prev_result.md", "k4_engine.py", "k4_quad.py"]:
    p = os.path.join(d, f)
    print("="*80)
    print(f, os.path.exists(p))
    if os.path.exists(p):
        with open(p) as fh:
            print(fh.read())
