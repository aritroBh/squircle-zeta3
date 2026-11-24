import json
base = "<project>"
with open(base + "/research/08_theory/SRF1_radii_3d/EV_coeffs_K30.json") as f:
    d = json.load(f)
print(list(d.keys()))
for k in d:
    print(k, len(d[k]))
print("V (first entries):")
for i, v in enumerate(d.get("V", [])[:8]):
    print(i, v[:40])
