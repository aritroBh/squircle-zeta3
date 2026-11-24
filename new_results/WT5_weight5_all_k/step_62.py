import json
W='<project>/research/08_theory/WT5_weight5_all_k/'
d1=json.load(open(W+'e3_lll_1100.json')); d2=json.load(open(W+'e3_lll_1300.json'))
def canon(v): return tuple(v if v[0]>0 else [-x for x in v])
r1=[x for x in d1['relations'] if x['vec'][0]!=0]
r2=[x for x in d2['relations'] if x['vec'][0]!=0]
c1={canon(x['vec']):x['residual'] for x in r1}
c2={canon(x['vec']):x['residual'] for x in r2}
common=set(c1)&set(c2)
print('r1',len(r1),'r2',len(r2),'identical',len(common))
for v in sorted(common):
    nz=[(i,c) for i,c in enumerate(v) if c]
    print('nonzero count', len(nz), 'resid', c1[v], c2[v])
    print('  ', nz[:20])
