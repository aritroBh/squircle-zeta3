import json
W='<project>/research/08_theory/WT5_weight5_all_k/'
d=json.load(open(W+'cobracket_e.json'))
print('W5pp_target:', d['W5pp_target'])
for nm in ['E1','E2','E3']:
    e=d[nm]
    print('='*10, nm, ' subset', e['subset'])
    print(' target:', e['target'])
    for t in e['correction']:
        print('   ', t['coefficient'], '* Li32(', t['x'], ',', t['y'], ')')
