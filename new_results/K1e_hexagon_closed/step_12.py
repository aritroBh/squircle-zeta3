
import os, json
base='research/08_theory'
k1c=os.path.join(base,'K1c_engine_predict')
print(os.listdir(k1c))
geo=json.load(open(os.path.join(k1c,'geometry.json')))
print(json.dumps(geo, indent=1)[:3500])
