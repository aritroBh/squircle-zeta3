base = "<project>/research/"
with open(base + "10_frontier/W1_corner_family/B_taylor_pi4/mu3_symbolic.py") as f:
    t = f.read()
i = t.find("def eval_formula")
print(t[i:i+2600])
