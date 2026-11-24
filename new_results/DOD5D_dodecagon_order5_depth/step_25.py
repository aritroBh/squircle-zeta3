import importlib, sys
sys.path.insert(0, ".")
import dod_engine
importlib.reload(dod_engine)
res170, diff170 = dod_engine.run(170, "p5_run170.txt")
import pickle
pickle.dump({"res": [str(x) for x in res170], "diff": str(diff170)}, open("p5_170.pkl","wb"))
