import json
from flint import arb, acb, ctx
ctx.prec = 600
p='<project>/research/08_theory/M5b_w6_closed_form/w6_600_partial.json'
d=json.load(open(p))
v=d['left']
x=acb(arb(v[0], arb(v[1])), arb(v[2], arb(v[3])))
print('parsed:', x.real.rad().str(8), x.is_finite())
print(x.real.mid().str(30))
