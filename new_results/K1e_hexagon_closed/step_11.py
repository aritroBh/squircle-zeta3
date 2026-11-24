
import os, json
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
geo=json.load(open(os.path.join(k1d,'geometry.json')))
print(json.dumps(geo, indent=1)[:4000])
