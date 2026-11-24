import os, json
base = "<project>/research"
cor4 = os.path.join(base, "08_theory/COR4_M3alpha_classical")
with open(os.path.join(cor4, "F_symbol.json")) as f:
    Fsym = json.load(f)
print(type(Fsym), len(Fsym))
# print a few entries
for i, (k, v) in enumerate(Fsym.items()):
    print(k, "->", v)
    if i > 5: break
