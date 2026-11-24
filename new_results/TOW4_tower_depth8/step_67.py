import pickle, os
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
K7 = pickle.load(open(f'{wts}/kern_7.pkl','rb'))
print("keys:", list(K7.keys()))
k50 = K7[(5,0)]
print("type (5,0):", type(k50), list(k50.keys()) if isinstance(k50, dict) else '')
for kk, vv in k50.items():
    print(" ", kk, type(vv), (list(vv.keys())[:8] if isinstance(vv, dict) else vv))
