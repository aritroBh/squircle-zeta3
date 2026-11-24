import os, json
base='<project>/research/08_theory/M4f_w5_extras'
d=json.load(open(os.path.join(base,'R5_min3.json')))
terms=d['terms']
li5=[k for k in terms if k.startswith('Li5')]
print('Li5 args in R5_min3:', li5)
print('other keys sample:', [k for k in terms if not k.startswith('Li5')][-15:])
print('N terms:', len(terms))
for k in li5: print(k, terms[k])
print('z5:', terms.get('z5'), 'Lchi:', terms.get('Lchi8_5'), [k for k in terms if 'chi' in k.lower() or 'L(' in k])
print()
print(open(os.path.join(base,'pool.json')).read()[:2500])
