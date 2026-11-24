import os, re
base = "<project>/research/08_theory/"
D1 = base+"DOD1_dodecagon/"
hits = []
for f in sorted(os.listdir(D1)):
    if f.startswith("step_") and f.endswith(".py"):
        txt = open(D1+f).read()
        if "'Lam'" in txt or '"Lam"' in txt or "Lam[" in txt or "Lam =" in txt or "def lam" in txt.lower():
            hits.append(f)
print(hits)
# find where atoms['Lam'] is created
for f in hits:
    txt = open(D1+f).read()
    if "atoms['Lam']" in txt or "'Lam':" in txt or "Lam={}" in txt or "Lam = {}" in txt:
        print("=====", f)
        print(txt[:4000])
        break
