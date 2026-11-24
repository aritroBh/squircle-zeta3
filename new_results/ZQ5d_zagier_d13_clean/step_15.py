import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
r8 = json.load(open(c+'/planner/d13_row8.json'))
print('d13_row8 keys:', list(r8.keys()))
for k in r8:
    v = r8[k]
    if isinstance(v, dict):
        print(k, ': dict with', len(v), 'entries; sample:', dict(list(v.items())[:3]))
    elif isinstance(v, list):
        print(k, ': list len', len(v), v[:5])
    else:
        print(k, ':', str(v)[:300])
