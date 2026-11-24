base = "<project>/research/"
import json, pickle
F4 = json.load(open(base + "08_theory/COR5_corner_classicality/F4_relation.json"))
print(json.dumps(F4, indent=1)[:3000])
