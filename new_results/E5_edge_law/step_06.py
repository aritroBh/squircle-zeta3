import os, json
d17='research/08_theory/D17_zeta5_polygons'
print(os.listdir(d17))
uw=json.load(open(os.path.join(d17,'universal_weights.json')))
print(json.dumps(uw, indent=1)[:3000])
