import os, json
print(os.getcwd())
print(os.listdir('.'))
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
with open(os.path.join(theory, 'genfun.json')) as f:
    genfun = json.load(f)
print(json.dumps(genfun['3|3,0'], indent=1)[:3000])
