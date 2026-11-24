import json, os
d = '<project>/research/08_theory/COR6_true_point_symbols'
print(os.listdir(d))
tbl = json.load(open(os.path.join(d, 'cor6_table.json')))
print(json.dumps(tbl, indent=1)[:6000])
