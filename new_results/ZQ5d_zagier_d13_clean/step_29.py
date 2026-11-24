import json
d8 = json.load(open('planner/d13_row8.json'))
prods = d8['products']
# check for same-multiset-different-order keys
from collections import defaultdict
canon = defaultdict(list)
for k in prods:
    canon[tuple(sorted(k.split('*')))].append(k)
dups = {c:ks for c,ks in canon.items() if len(ks)>1}
print('product keys:', len(prods), 'canonical:', len(canon), 'dup-canonical groups:', len(dups))
for c,ks in list(dups.items())[:10]:
    print(' ', ks, '->', [prods[k] for k in ks])
# also stats of row8 json
li5 = d8['li5']
from fractions import Fraction as Ff
allc = [Ff(v) for v in list(li5.values())+list(prods.values())+[d8['z5_coef']]]
import math
print('row8: nLi5', len(li5), 'nprod', len(prods), 'max|num|', max(abs(x.numerator) for x in allc),
      'lcm', math.lcm(*[x.denominator for x in allc]))
print('T_coef', d8['T_coef'], '=', Ff(d8['T_coef'])/13**5, '*13^5')
