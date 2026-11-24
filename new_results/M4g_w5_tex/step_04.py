import os, json
base='../M4e_w5_lchi'
print(os.path.exists(base+'/planner/R5_nice_reduced.json'))
print(os.path.exists('../../../paper/ieee/main.tex'))
d=json.load(open(base+'/planner/R5_nice_reduced.json'))
print(list(d.keys()))
terms=d['terms']
print(type(terms), len(terms))
for k,v in terms.items():
    print(repr(k), v)
