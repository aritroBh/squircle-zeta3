import json
fs = json.load(open('<project>/research/08_theory/COR4_M3alpha_classical/F_symbol.json'))
print(type(fs))
if isinstance(fs, dict):
    print(list(fs.keys())[:10])
    s=json.dumps(fs)
    print(s[:1500])
elif isinstance(fs, list):
    print(len(fs))
    print(json.dumps(fs[:3], indent=1)[:1500])
