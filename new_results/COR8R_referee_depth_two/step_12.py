import json, os
base = '<project>'
sa = json.load(open(os.path.join(base,'research/08_theory/COR4_M3alpha_classical/sunit_args.json')))
print(type(sa))
if isinstance(sa, dict):
    print(list(sa.keys()))
elif isinstance(sa, list):
    print(len(sa))
    print(json.dumps(sa[:5], indent=1))
