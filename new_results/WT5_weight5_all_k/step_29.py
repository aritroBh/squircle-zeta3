base='<project>'
with open(base+'/research/08_theory/WT4_weight4_all_k/progress.log') as f:
    print(f.read())
print('---- quad5_100 ----')
import json
print(json.load(open(base+'/research/08_theory/M5b_w6_closed_form/quad5_100.json')))
print('---- Lchi8_5.json digits ----')
v=json.load(open(base+'/research/08_theory/M5b_w6_closed_form/Lchi8_5.json'))
print({k:(len(str(vv)) if isinstance(vv,str) else vv) for k,vv in v.items()})
v2=json.load(open(base+'/research/08_theory/M5b_w6_closed_form/Lchi8_5_1800_planner.json'))
print({k:(len(str(vv)) if isinstance(vv,str) else vv) for k,vv in v2.items()})
