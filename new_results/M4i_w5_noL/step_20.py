import json
from fractions import Fraction as F
dst = '<project>/research/08_theory/M4i_w5_noL'
d = json.loads(open(dst+'/L_relation.json').read())
terms = d['terms']
print('n terms:', len(terms))
for nm,(p,q) in terms.items():
    print(f'  {nm:28s} [{p}, {q}]')
# true 10 largest
allc = []
for nm,(p,q) in terms.items():
    mx = max(abs(F(p).numerator), abs(F(q).numerator))
    allc.append((mx, nm, p, q))
allc.sort(reverse=True)
print('\n10 largest:')
for mx,nm,p,q in allc[:10]:
    print(f'  |num|={mx}  {nm}: [{p}, {q}]')
from math import lcm
dens = [F(x).denominator for pq in terms.values() for x in pq]
nums = [abs(F(x).numerator) for pq in terms.values() for x in pq]
print('max|num| =', max(nums), ' lcm =', lcm(*dens), ' terms =', len(terms))
