import json
m5b='<project>/research/08_theory/M5b_w6_closed_form'
d=json.load(open(m5b+'/w5pp_4800.json'))
print({k:v for k,v in d.items() if 'value' not in k})
p=json.load(open(m5b+'/w5pp_4800_partial.json'))
for k in ['left','mid','right']:
    print(k, 'rad:', p[k][1])
log=open(m5b+'/progress.log').read().splitlines()
print('\n'.join(log[-12:]))
