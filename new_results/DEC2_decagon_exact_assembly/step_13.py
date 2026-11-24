p = "<project>/paper/ieee/main.tex"
with open(p) as f:
    tex = f.read()
for needle in ["1099", "1991", "S_5", "L_5", "81/64", "tabular"]:
    idxs = []
    i = tex.find(needle)
    while i != -1:
        idxs.append(i)
        i = tex.find(needle, i+1)
    print(needle, idxs[:10])
