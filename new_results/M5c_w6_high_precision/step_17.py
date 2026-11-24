import json, shutil, os
base='<project>/research/08_theory/M5b_w6_closed_form'
work='.'
for s in ['omega5.json','quad5_100.json','W6_400.txt','w6_8000.json']:
    shutil.copy(os.path.join(base,s), os.path.join(work,s))
pool=json.load(open('pool.json'))
print(json.dumps(pool, indent=1)[:2000])
n6=json.load(open('n6_solution.json'))
print(json.dumps(n6, indent=1))
