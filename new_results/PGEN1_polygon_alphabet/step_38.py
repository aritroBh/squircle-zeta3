import os, pickle
base = "<project>"
O6 = os.path.join(base,"research/08_theory/OCT6_octagon_closed")
t = pickle.load(open(os.path.join(O6,"table_Stot.pkl"),"rb"))
print(type(t))
for k,v in (t.items() if isinstance(t,dict) else enumerate(t)):
    print(repr(k)[:60], "=>", repr(v)[:120])
