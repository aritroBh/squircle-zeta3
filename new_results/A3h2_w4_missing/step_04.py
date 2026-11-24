import json
base='<project>'
# Print section 7 of A3f4 result.md
txt=open(base+'/research/08_theory/A3f4_glanois_basis/result.md').read()
i=txt.find('## 7')
print(txt[i:i+2500] if i>=0 else 'no section 7 header')
print('====T3_quotient_nf.json====')
d=json.load(open(base+'/research/08_theory/A3h_weight4_target/T3_quotient_nf.json'))
print(type(d), list(d.keys()) if isinstance(d,dict) else len(d))
