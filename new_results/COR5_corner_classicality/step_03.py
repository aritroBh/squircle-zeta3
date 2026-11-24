import os
base = "<project>/research"
cor4 = os.path.join(base, "08_theory/COR4_M3alpha_classical")
with open(os.path.join(cor4, "result.md")) as f:
    txt = f.read()
print(txt[:14000])
