import json
W='<project>/research/08_theory/WT5_weight5_all_k/'
for t in ['e1','e2','e3']:
    v=open(W+f'{t}_1500.txt').read().strip()
    print(t, v[:70], '... len', len(v))
d=json.load(open(W+'check_tsv_vs_quad.json'))
print(d)
