import os, json
target = '<project>/research/08_theory/COR8B_depth_two_two_terms'
with open(os.path.join(target, 'twovar_symbols.json')) as f:
    tws = json.load(f)
print(json.dumps(tws, indent=1)[:6000])
