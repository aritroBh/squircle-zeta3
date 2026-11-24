import json
WORK='<project>/research/08_theory/M4f_w5_extras'
T='<project>/research/08_theory'
d = json.loads(open(f'{WORK}/R5_best.json').read())['terms']
prods = [k for k in d if '*' in k]
print(len(prods), 'product terms in R5_best:')
for k in sorted(prods): print('  ', k, d[k])
print()
w4 = json.loads(open(f'{T}/M4a_decomposition/w4pp_closed_form.json').read())
print(json.dumps(w4, indent=1)[:800])
