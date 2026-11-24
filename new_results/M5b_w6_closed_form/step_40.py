import json, mpmath as mp
mp.mp.dps=200
d=json.load(open('<project>/research/08_theory/M5b_w6_closed_form/w6_600.json'))
print(d['value'][:80])
ref=open('<project>/research/08_theory/M5b_w6_closed_form/W6_400.txt').read().strip()
print(ref[:80])
v=mp.mpf(d['value']); r=mp.mpf(ref)
print('diff', mp.nstr(abs(v-r),5))
print('radius in json:', d['radius'], 'tail:', d['tail_bound'])
# sum of checkpoint pieces
p=json.load(open('<project>/research/08_theory/M5b_w6_closed_form/w6_600_partial.json'))
tot=mp.fsum([mp.mpf(p[k][0]) for k in ['left','mid','right']])
print('sum of ckpt pieces:', mp.nstr(tot,40))
print('json value        :', mp.nstr(v,40))
print('diff json - ckptsum:', mp.nstr(abs(v-tot),5))
