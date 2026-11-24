import json
M4E = '<project>/research/08_theory/M4e_w5_lchi'
d = json.load(open(M4E+'/lchi_search_1200.json'))
cn = d['colnames']
h0 = d['history'][0]
print("round0 deps:")
for dep in h0['deps']:
    vec = dep['vec']
    sup = [(cn[j], c) for j, c in enumerate(vec) if c]
    print(" maxcoef", dep['maxcoef'], "resid", dep['residual'])
    print("   ", sup)
