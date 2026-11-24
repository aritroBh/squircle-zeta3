import json, shutil, os
from collections import Counter
base = '<project>/research/08_theory'
shutil.copy(f'{base}/ZQ5b_zagier_w5_larger_S/d17_final.json', os.getcwd())
J = json.load(open('d17_final.json'))
prods = J['products']
print("D=17 products:", len(prods))
types = Counter()
for k in prods:
    parts = k.split('*')
    npi2 = parts.count('pi2'); npi4 = parts.count('pi4')
    nz3 = sum(1 for p in parts if p=='z3'); nt3 = sum(1 for p in parts if p=='T3')
    logs = [p for p in parts if p not in ('pi2','pi4','z3','T3')]
    key = (npi2, npi4, nz3, nt3, len(logs))
    types[key]+=1
print("(pi2,pi4,z3,T3,nlogs) -> count:")
for k,v in sorted(types.items()):
    print("  ", k, v)
print("\nletters:", J['letters'])
print("z5_coef:", J['z5_coef'], "T_coef:", J['T_coef'])
for i,k in enumerate(list(prods)[:15]):
    print("  ", k, "=", prods[k])
