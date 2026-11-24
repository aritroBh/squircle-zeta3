base = '<project>/research/'
import json
tab = json.load(open(base+'08_theory/COR6_true_point_symbols/cor6_table.json'))
print(json.dumps(tab, indent=1))
