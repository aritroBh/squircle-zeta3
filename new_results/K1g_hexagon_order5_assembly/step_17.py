import os, json
base = '<project>/research/08_theory'
k1c = os.path.join(base, 'K1c_engine_predict')
with open(os.path.join(k1c, 'geometry.json')) as f:
    geo = json.load(f)
import pprint
print("=== geometry['0'] (b=0 hexagon) ===")
pprint.pprint(geo['0'])
print()
print("=== geometry['1'] ===")
pprint.pprint(geo['1'])
