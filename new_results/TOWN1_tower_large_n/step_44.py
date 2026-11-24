import os, json
base = '<project>'
p = os.path.join(base, 'research/08_theory/WT4_weight4_all_k')
with open(os.path.join(p, 'kernels_deg3.json')) as f:
    kd3 = json.load(f)
print(type(kd3), list(kd3.keys())[:10] if isinstance(kd3, dict) else len(kd3))
print(json.dumps(kd3, indent=1)[:3500])
