import json
base='research/'
d = json.load(open(base+'10_frontier/H1_octahedron_end/routeA_n2_K12.json'))
v = d['2']
print(v['d'][:5])
src = open(base+'08_theory/D11c_zeta_first_order/step_02.py').read()
print(src[:5500])
