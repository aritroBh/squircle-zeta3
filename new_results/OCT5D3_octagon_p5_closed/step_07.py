import json, pickle
base = "<project>/research/08_theory/"
with open(base + "M4e_w5_lchi/result.md") as f:
    txt = f.read()
print(txt[:12000])
