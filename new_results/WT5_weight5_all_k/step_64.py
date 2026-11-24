import json
W='<project>/research/08_theory/WT5_weight5_all_k/'
d=json.load(open(W+'e3_closed_form_rel0.json'))
for nm,a,b in zip(d['names'],d['ca'],d['cb']):
    print(f'{a} + ({b})*sqrt2   * {nm}')
print('aL,bL:',d['aL'],d['bL'])
