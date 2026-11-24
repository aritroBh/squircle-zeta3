import json
base = '<project>/research/08_theory'
s1 = f'{base}/SUB1_subleading_table'
tab = json.load(open(f'{s1}/subtop_table.json'))
print(type(tab), list(tab.keys())[:10])
import pprint
e3 = tab.get('3') or tab.get(3)
pprint.pprint(str(e3)[:1500])
