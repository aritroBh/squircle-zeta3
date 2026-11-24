import json, os
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
with open(os.path.join(theory, 'genfun.json')) as f:
    genfun = json.load(f)
print(type(genfun))
if isinstance(genfun, dict):
    print(list(genfun.keys())[:20])
elif isinstance(genfun, list):
    print(len(genfun))
    print(json.dumps(genfun[0], indent=1)[:2000])
