import json
from fractions import Fraction as F
from math import lcm
dst = '<project>/research/08_theory/M4i_w5_noL'
def stats(path):
    terms = json.loads(open(dst+'/'+path).read())['terms']
    nums = [abs(F(x).numerator) for pq in terms.values() for x in pq]
    dens = [F(x).denominator for pq in terms.values() for x in pq]
    li5 = [nm for nm in terms if nm.startswith('Li5')]
    prods = [nm for nm in terms if '*' in nm]
    return dict(n=len(terms), mx=max(nums), lcm=lcm(*dens), li5=len(li5), prods=len(prods),
                hasL='Lchi8_5' in terms, z5='z5' in terms)
for p in ['R5_best.json','R5_noL.json','L_relation.json']:
    print(p, stats(p))
print()
t = json.loads(open(dst+'/L_relation.json').read())['terms']
print('L-relation, full term list (name: [p, q] with L = sum (p+q*sqrt2)*atom):')
for nm,(p,q) in t.items():
    print(f'  {nm:28s} [{p}, {q}]')
print()
print(open(dst+'/progress.log').read())
import os
steps = sorted(f for f in os.listdir(dst) if f.startswith('step_') and f.endswith('.py'))
print('step files:', steps)
