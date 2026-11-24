import os, json
base = '<project>/research/08_theory'
k1c = os.path.join(base, 'K1c_engine_predict')
print("=== K1c contents ===")
for d in sorted(os.listdir(k1c)):
    print(d)
print()
gj = os.path.join(k1c, 'geometry.json')
if os.path.exists(gj):
    with open(gj) as f:
        geo = json.load(f)
    print("=== geometry.json keys ===")
    print(list(geo.keys()))
