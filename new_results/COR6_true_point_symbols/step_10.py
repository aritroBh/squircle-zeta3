base = "<project>/research/"
p = base + "08_theory/COR4_M3alpha_classical/"
with open(p + "cor4_sym.py") as f:
    src = f.read()
print(src[6000:])
