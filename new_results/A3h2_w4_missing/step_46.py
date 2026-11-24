import sys, pickle, json
from math import gcd
from functools import reduce
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
d=pickle.load(open(WORK+'/cands29_reduced.pickle','rb'))
fw=d['freewords']; red=d['cands']
# make primitive + normalize sign (leading coeff positive)
def normalize(row):
    g=reduce(gcd,[abs(x) for x in row if x])
    row=[x//g for x in row]
    for x in row:
        if x: 
            if x<0: row=[-y for y in row]
            break
    return row
red=[normalize(r) for r in red]
heights=[max(abs(x) for x in r) for r in red]
print('heights:',sorted(heights))
print('gcd check ok; total nonzero coeffs:',sum(sum(1 for x in r if x) for r in red))
# relation 0 explicitly
print('\n=== relation N0 (form-level, support 6) ===')
for i,x in enumerate(red[0]):
    if x: print(f'  {x:+d} * {fw[i]}')
# save JSON: list of relations, each {word: coeff}
out={'freewords':[list(w) for w in fw],
     'relations':[{ '.'.join(str(a) for a in fw[i]): c for i,c in enumerate(r) if c} for r in red],
     'meta':'29 missing weight-4 level-8 relations; verified 680 digits; heights %d..%d'%(min(heights),max(heights))}
json.dump(out,open(WORK+'/w4_missing_relations.json','w'))
print('\nsaved w4_missing_relations.json')
# also verify relation 0 at all 4 embeddings with stored 350-digit values
import mpmath as mp
mp.mp.dps=340
out2=pickle.load(open(WORK+'/vals110.pickle','rb'))
V={k:[mp.mpc(mp.mpf(out2['vals'][k][w][0]),mp.mpf(out2['vals'][k][w][1])) for w in fw] for k in (1,3,5,7)}
for k in (1,3,5,7):
    s=mp.mpc(0)
    for i,c in enumerate(red[0]):
        if c: s+=c*V[k][i]
    print(f'N0 at k={k}: {mp.nstr(abs(s),3)}')
