import json
p = '<project>/research/08_theory/M4e_w5_lchi/transcript.json'
d = json.load(open(p))
print(type(d))
if isinstance(d, dict):
    print(list(d.keys()))
elif isinstance(d, list):
    print(len(d))
    print(json.dumps(d[:2], indent=2)[:3000])
