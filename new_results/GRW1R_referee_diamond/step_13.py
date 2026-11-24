import json, os
p = '../D11b2_stokes_jump/step2_num450_N180.json'
print(os.path.exists(p))
d = json.load(open(p))
print(type(d), list(d.keys())[:10] if isinstance(d,dict) else len(d))
f450 = d['f450']
print(len(f450), type(f450[0]), f450[0][:60], f450[1][:60])
print(d.get('meta', {}) if isinstance(d,dict) else '')
print({k:v for k,v in d.items() if k!='f450'} if isinstance(d,dict) else '')
