
import os, pickle, json
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
WTS = os.path.join(T, 'WTS_weight_k_structure/theory')
K5 = pickle.load(open(os.path.join(WTS, 'kern_5.pkl'), 'rb'))
print("kern_5.pkl type:", type(K5), "keys:", list(K5.keys()) if hasattr(K5, 'keys') else None)
for k, v in K5.items():
    print(f"\n--- key {k} ---")
    print(type(v))
    print(repr(v)[:1500])
