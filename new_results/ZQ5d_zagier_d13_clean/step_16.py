import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
d17 = json.load(open(c+'/d17_final.json'))
print('d17_final keys:', list(d17.keys()))
for k in d17:
    v=d17[k]
    if isinstance(v,dict):
        print(k,': dict',len(v),'sample:',dict(list(v.items())[:2]))
    else:
        print(k,':',str(v)[:200])
# check ratio row8
from fractions import Fraction as Ff
print(Ff(-13216098208787298968, 35594794969976))
import math
print('gcd check:', math.gcd(-13216098208787298968, 35594794969976))
print(-13216098208787298968 // 371293, 371293*35594794969976)
