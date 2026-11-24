import json
p='<project>/research/08_theory/M5b_w6_closed_form/w6_600_partial.json'
d=json.load(open(p))
print(d['left'])
from flint import arb, acb
v=d['left']
try:
    x=acb(arb(v[0], arb(v[1])), arb(v[2], arb(v[3])))
    print('parsed:', x.str(10), x.is_finite())
except Exception as e:
    print('parse error:', e)
    print('arb rad test:', arb(v[1]))
