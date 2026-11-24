
base = "<project>/research/08_theory"
import json
d = json.load(open(base + "/COR5_corner_classicality/F4_relation.json"))
print(json.dumps(d, indent=1)[:4000])
