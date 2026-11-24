import json
base = "<project>"
with open(base + "/research/08_theory/SRF1_radii_3d/SRF1_n3_K60.json") as f:
    d60 = json.load(f)
print(type(d60), list(d60.keys()) if isinstance(d60, dict) else len(d60))
s = json.dumps(d60)
print(s[:800])
