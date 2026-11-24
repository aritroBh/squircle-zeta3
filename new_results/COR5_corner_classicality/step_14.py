import os
base = "<project>/research"
g1 = os.path.join(base, "10_frontier/G1_zeta3_any_angle")
with open(os.path.join(g1, "result.md")) as f:
    txt = f.read()
# find section B3
i = txt.find("B3")
print(txt[i-2000:i+9000])
