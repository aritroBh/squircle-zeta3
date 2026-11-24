
base = "<project>/research/08_theory"
import json
d = json.load(open(base + "/M4c_w5_closed_form/cobracket.json"))
print("W5 cobracket:", d["W5"])
print()
for k, v in d["depth_two"].items():
    print(k, "->", v)
print()
print(json.dumps(d["correction"], indent=1))
