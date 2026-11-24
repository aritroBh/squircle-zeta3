base='<project>'
import json
om=json.load(open(base+'/research/08_theory/M5b_w6_closed_form/omega5.json'))
print(list(om.keys()))
for k in om:
    print(k, om[k])
