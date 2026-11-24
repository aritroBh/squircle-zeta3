base = "<project>/research/08_theory/"
with open(base+"OCT9b_octagon_debug/result.md") as f:
    txt = f.read()
print(len(txt))
print(txt[:9000])
