import os, json
base = '<project>'
p = os.path.join(base, 'research/08_theory/WT4_weight4_all_k')
for fn in ['closed_forms.json']:
    with open(os.path.join(p, fn)) as f:
        data = json.load(f)
    print("="*30, fn)
    print(json.dumps(data, indent=1)[:4000])
