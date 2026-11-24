
import os, json
base='research/08_theory'
k1c=os.path.join(base,'K1c_engine_predict')
t3=json.load(open(os.path.join(k1c,'t3_b0_150.json')))
print(t3.keys())
print(json.dumps({k:v for k,v in t3.items() if k!='arcs'}, indent=1)[:2000])
arc0=t3['arcs'][0]
print(arc0.keys())
for k,v in arc0.items():
    print(k, str(v)[:200])
