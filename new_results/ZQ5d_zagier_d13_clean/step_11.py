import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
k5 = json.load(open(c+'/out/d13/kernel5_cap.json'))
print('kernel5_cap keys:', list(k5.keys()) if isinstance(k5,dict) else type(k5))
if isinstance(k5, dict):
    for k,v in k5.items():
        print(k, type(v), (len(v) if hasattr(v,'__len__') else v))
