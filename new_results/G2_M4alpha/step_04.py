import os
print(os.getcwd())
print(sorted(os.listdir('.')))
base = "<project>/research/"
with open(base + "10_frontier/G1_zeta3_any_angle/result.md") as f:
    txt = f.read()
print(len(txt))
print(txt[12000:20000])
