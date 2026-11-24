import pickle, importlib, sys
sys.path.insert(0, ".")
P3_630 = pickle.load(open("<project>/research/08_theory/DOD1_dodecagon/P3_630.pkl","rb"))
print(P3_630.keys())
for k,v in P3_630.items():
    print(k, str(v)[:70])
import dod_engine
importlib.reload(dod_engine)
res120, diff120 = dod_engine.run(120, "p5_run120.txt")
