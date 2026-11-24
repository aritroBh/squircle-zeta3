base = "<project>/research/"
p = base + "08_theory/COR4_M3alpha_classical/"
import os
print(os.listdir(p))
with open(p + "cor4_sym.py") as f:
    print(f.read()[:6000])
