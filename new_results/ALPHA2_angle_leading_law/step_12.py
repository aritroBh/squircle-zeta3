import pickle
with open('beta_squircle.pkl','rb') as f: bs = pickle.load(f)
from fractions import Fraction
# print squircle table organized by d
ds = sorted(set(k[0] for k in bs))
for d in ds:
    print(f"--- d={d} ---")
    keys = sorted(set((k[1],k[2]) for k in bs if k[0]==d))
    for kk in keys:
        print(kk, bs[(d,)+kk])
