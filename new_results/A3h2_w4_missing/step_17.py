import json, pickle
base='<project>'
# Rogers five-term
d=json.load(open(base+'/research/08_theory/A3f_level8_relations/weight2_new_relation.json'))
print('weight2_new_relation.json type:', type(d))
if isinstance(d,dict): print(list(d.keys()))
print(str(d)[:1500])
