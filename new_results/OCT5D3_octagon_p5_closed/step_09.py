import json
base = "<project>/research/08_theory/"
r5 = json.load(open(base+"M4f_w5_extras/R5_best.json"))
terms = r5["terms"]
print("number of terms:", len(terms))
for k,v in terms.items():
    print(k, v)
