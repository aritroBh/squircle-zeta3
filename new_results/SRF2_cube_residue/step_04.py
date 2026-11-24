import json
base = "<project>"
with open(base + "/research/08_theory/SRF1_radii_3d/EV_coeffs_K30.json") as f:
    d = json.load(f)
print(json.dumps(d, indent=1)[:6000])
