import json
M4D = '<project>/research/08_theory/M4d_w5_nice'
d = json.load(open(M4D+'/nice_found_1200.json'))
print(list(d.keys()))
print("digits:", d['digits'])
print("extras:", d['extras'])
print("vec:", d['vec'])
print("residual:", d['residual'])
print("num K rows:", len(d['K']))
