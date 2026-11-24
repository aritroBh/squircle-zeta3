import json
r = json.load(open('../COR9_golden_corner_closed_forms/gol1_W4_relation.json'))
print(json.dumps(r, indent=1)[:3000])
