import json
M4E = '<project>/research/08_theory/M4e_w5_lchi'
d = json.load(open(M4E+'/lchi_sat_search_1200.json'))
cn = d['colnames']
for h in d['history']:
    print("round", h['round'], "kept", len(h['kept']), "genuine", len(h['genuine']), "deps", len(h['deps']))
    for dep in h['deps']:
        sup = [(cn[j], c) for j, c in enumerate(dep['vec']) if c]
        print("   dep maxcoef", dep['maxcoef'], "resid", dep['residual'], sup)
