import json
W='<project>/research/08_theory/WT5_weight5_all_k/'
for t in ['e1','e2','e3','w5pp']:
    d=json.load(open(W+f'{t}_600.json'))
    print(t, d['value'][:60], ' radius', d['radius'])
c=json.load(open(W+'check_tsv_vs_quad.json'))
for k,v in c.items():
    print(k, v['err_vs_quad'], ' G=', v['G_value'][:50])
