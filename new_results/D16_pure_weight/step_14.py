import json, os
base='research/'
d = json.load(open(base+'10_frontier/H1_octahedron_end/routeA_n2_K12.json'))
print(type(d), list(d.keys()))
v = d['2'] if '2' in d else d[list(d.keys())[0]]
print(v.keys())
print(json.dumps(v, indent=1)[:1500])
