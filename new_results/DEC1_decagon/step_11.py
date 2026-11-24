import json
d='research/08_theory/D15b_polygon_law/'
v=json.load(open(d+'vals_m2_8_200.json'))
print(type(v), list(v.keys()) if isinstance(v,dict) else len(v))
entry=v['5']
print(type(entry), len(entry))
for i,e in enumerate(entry):
    print(i, str(e)[:120])
