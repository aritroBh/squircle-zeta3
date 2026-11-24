import json
base='research/08_theory/M4e_w5_lchi'
d=json.load(open(base+'/planner/R5_nice_reduced.json'))
print(list(d.keys()))
terms=d['terms']
print(type(terms), len(terms))
for k,v in terms.items():
    print(repr(k), v)
