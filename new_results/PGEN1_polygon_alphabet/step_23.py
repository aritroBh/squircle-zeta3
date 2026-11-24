import os, pickle, sys
base = "<project>"
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
print(type(sf), list(sf.keys()) if isinstance(sf, dict) else len(sf))
import pprint
if isinstance(sf, dict):
    for k,v in sf.items():
        print("KEY:", k, type(v))
        s = repr(v)
        print(s[:1500])
        print("---")
