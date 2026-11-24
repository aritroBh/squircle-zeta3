base = '<project>/research/'
src = open(base+'08_theory/COR5_corner_classicality/result.md').read()
print(src[4300:6100])
print('='*80)
import json
print(json.dumps(json.load(open(base+'08_theory/COR5_corner_classicality/F4_relation.json')), indent=1))
