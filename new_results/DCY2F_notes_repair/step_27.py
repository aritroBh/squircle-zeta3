import os, re
base = "<project>/research/08_theory/DCY2_decay_first_correction"
pats = ["2*bR*y**3", "2*b_R*y**3", "h_1+", "h1 +", "3/2 y", "1.5y", "E_tailmom", "2*Et", "1409", "0.00141"]
for fn in sorted(os.listdir(base)):
    if not fn.endswith((".py", ".md")): continue
    with open(os.path.join(base, fn)) as f:
        txt = f.read()
    hits = []
    for p in pats:
        if p in txt: hits.append(f"{p}({txt.count(p)})")
    if hits:
        print(f"{fn}: {', '.join(hits)}")
