import json
from fractions import Fraction as F
pool=json.load(open('pool.json'))['pool']
for i,e in enumerate(pool):
    print(i, e['p'], e['q'], 'h=',e['h'])
