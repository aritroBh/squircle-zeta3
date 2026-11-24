base = "<project>/research/08_theory/"
with open(base+"OCT9c_octagon_close/result.md") as f:
    txt = f.read()
print(len(txt))
print(txt[6000:14000])
