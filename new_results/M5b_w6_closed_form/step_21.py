import os, json
base='<project>/research/08_theory/M4f_w5_extras'
d=json.load(open(os.path.join(base,'R5_min3.json')))
print(type(d), list(d.keys()) if isinstance(d,dict) else len(d))
print(json.dumps(d, indent=1)[:3000])
