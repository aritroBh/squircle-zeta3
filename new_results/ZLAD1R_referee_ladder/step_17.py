import json
base = '<project>'
d = base + '/research/08_theory/ZLAD1_zeta_ladder'
t8 = json.load(open(d + '/mu_exact_N150_K8.json'))
print('keys k:', sorted(t8.keys(), key=int))
for k in ['7','8']:
    print('k=%s n range:'%k, min(map(int,t8[k].keys())), max(map(int,t8[k].keys())))
print('mu_{8,1}:', t8['8']['1'])
print('mu_{7,2}:', t8['7']['2'])
t6 = json.load(open(d + '/mu_exact_N300_K6.json'))
print('K6 file: k keys', sorted(t6.keys(), key=int), 'max n', max(map(int, t6['6'].keys())))
