import json
base='../M4e_w5_lchi'
d=json.load(open(base+'/planner/R5_nice_reduced.json'))
print('source:', d['source'])
print('N5:', json.dumps(d['N5'], indent=1)[:2000])
c=json.load(open(base+'/cobracket.json'))
print(list(c.keys()))
print(json.dumps(c['correction'], indent=1))
