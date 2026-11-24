import os, pickle
base = "<project>"
DOD = os.path.join(base,"research/08_theory/DOD5_dodecagon_termwise")
t = pickle.load(open(os.path.join(DOD,"dod5_tallies.pkl"),"rb"))
print(type(t))
if isinstance(t, dict):
    for k,v in t.items():
        print("KEY", k, type(v), repr(v)[:2000])
