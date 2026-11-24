import json
WORK='<project>/research/08_theory/M4f_w5_extras'
d = json.loads(open(f'{WORK}/overcomplete_1200.json').read())
print("kept:", d['kept'])
print()
for dep in d['deps']:
    kept = dep['kept']
    sup = [(kept[j], c) for j,c in enumerate(dep['vec'][1:]) if c]  # vec[0] is target(=0)
    print(f"maxcoef={dep['maxcoef']} resid={dep['resid']}")
    print("   ", sup)
