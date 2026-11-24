import os, json
base = '<project>'
p = os.path.join(base, 'research/08_theory/WT3_weight3_all_k')
wk = json.load(open(os.path.join(p,'wt3_kernels.json')))
print(type(wk), list(wk.keys())[:5] if isinstance(wk, dict) else len(wk))
# inspect structure
import pprint
if isinstance(wk, dict):
    for kk in list(wk.keys())[:2]:
        print("KEY", kk)
        pprint.pprint(wk[kk])
