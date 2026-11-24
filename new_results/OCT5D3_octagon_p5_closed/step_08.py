import json
base = "<project>/research/08_theory/"
r5 = json.load(open(base+"M4f_w5_extras/R5_best.json"))
print(type(r5))
if isinstance(r5, dict):
    print(list(r5.keys()))
    for k,v in r5.items():
        print("---", k, type(v), str(v)[:800])
