p = '../PGEN1_polygon_alphabet/result.md'
with open(p) as f:
    txt = f.read()
print(len(txt))
print(txt[:8000])
