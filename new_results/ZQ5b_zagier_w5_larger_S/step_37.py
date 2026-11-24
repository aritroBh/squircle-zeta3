import os, json, pickle
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
P = '../ZQ5_zagier_weight5/planner/'
R = json.load(open(P+'d17_relation_planner.json'))
print("keys:", list(R.keys()))
for k,v in R.items():
    if isinstance(v,(str,int,float)): print(k, '=', v)
    elif isinstance(v,list): print(k, ': list len', len(v), ' first:', v[:2])
    elif isinstance(v,dict): print(k, ': dict keys', list(v.keys())[:10])
