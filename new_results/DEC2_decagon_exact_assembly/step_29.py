import os, re
base = "<project>/research/08_theory/OCT6_octagon_closed"
hits = []
for f in sorted(os.listdir(base)):
    if f.endswith(".py"):
        with open(os.path.join(base,f)) as fh:
            txt = fh.read()
        if "D3" in txt and ("quad" in txt or "integr" in txt.lower() or "def " in txt):
            hits.append(f)
print(hits)
