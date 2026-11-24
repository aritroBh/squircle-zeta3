import json
W='<project>/research/08_theory/WT5_weight5_all_k/'
for t in ['e1','e2','e3']:
    d=json.load(open(W+f'{t}_5800.json'))
    v=d['value_full']
    with open(W+f'{t}_1500.txt','w') as f:
        f.write(v[:v.find('.')+1501]+'\n')
    print(t, 'digits saved 1500 of', len(v), 'value head:', v[:50])
with open(W+'progress.log','a') as f:
    f.write('step_D hp values: e1/e2/e3 at 5400&5800 bits; agreements 3.007e-1595 / 2.773e-1595 / 3.020e-1595; e{1,2,3}_1500.txt saved\n')
print('logged')
