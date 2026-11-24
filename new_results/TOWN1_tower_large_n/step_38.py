import os, pickle
base = '<project>'
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
k4 = pickle.load(open(os.path.join(WTS,'kern_4.pkl'),'rb'))
print("kern_4 type:", type(k4), "keys:", list(k4.keys()) if isinstance(k4, dict) else len(k4))
for key in k4:
    print("key", key, "->", type(k4[key]))
    print("   ", k4[key])
    break
# Look at the degree-3 keys
for key in [(3,0),(2,1)]:
    if key in k4:
        print("deg3 key", key, ":", k4[key])
