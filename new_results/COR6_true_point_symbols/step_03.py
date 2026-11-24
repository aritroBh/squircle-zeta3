base = "<project>/research/"
with open(base + "08_theory/COR5_corner_classicality/result.md") as f:
    txt = f.read()
print(txt[12000:20000])
print("LEN", len(txt))
