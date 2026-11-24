import os, json, pickle
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
P = '../ZQ5_zagier_weight5/planner/'
R = json.load(open(P+'d17_relation_planner.json'))
print("letters:", R['letters'])
print("T_coef:", R['T_coef'], " z5_coef:", R['z5_coef'])
li5 = R['li5']
print("n li5:", len(li5))
items = list(li5.items())[:8]
for k,v in items: print("  ", k, "->", v)
print("...")
products = R['products']
print("n products:", len(products))
# the pkl with zeros
raw = pickle.load(open(P+'d17_rerun_1150.pkl','rb'))
print("\npkl keys:", list(raw.keys()))
for k,v in raw.items():
    if isinstance(v,list): print(k, "list len", len(v), "first:", str(v[0])[:120] if v else "")
    else: print(k, ":", str(v)[:150])
