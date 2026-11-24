import json
base = "<project>/research/08_theory/"
cb = json.load(open(base+"M4e_w5_lchi/cobracket.json"))
print(json.dumps(cb, indent=1)[:4000])
