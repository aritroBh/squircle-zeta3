import os
base = "<project>/research"
g1 = os.path.join(base, "10_frontier/G1_zeta3_any_angle")
print(os.listdir(g1))
with open(os.path.join(g1, "result.md")) as f:
    txt = f.read()
print("LEN:", len(txt))
