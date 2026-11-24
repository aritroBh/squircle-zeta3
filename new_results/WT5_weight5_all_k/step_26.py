base='<project>'
import json
pool=json.load(open(base+'/research/08_theory/M5b_w6_closed_form/pool.json'))
print(json.dumps(pool)[:2000])
print('---Lchi8---')
print(json.load(open(base+'/research/08_theory/M5b_w6_closed_form/Lchi8_5.json'))['value'][:80])
print('---progress tail---')
with open(base+'/research/08_theory/M5b_w6_closed_form/progress.log') as f:
    lines=f.read().splitlines()
for l in lines[-30:]: print(l)
