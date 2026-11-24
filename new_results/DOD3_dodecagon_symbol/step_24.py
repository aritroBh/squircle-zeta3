import sys, os, pickle
sys.path.insert(0, os.getcwd())
base = os.path.dirname(os.getcwd())
d2 = os.path.join(base,'DOD2_dodecagon_finish')
with open(os.path.join(d2,'R256_sym.pkl'),'rb') as f:
    R = pickle.load(f)
print(type(R))
if isinstance(R, dict):
    for k,v in R.items():
        print('---', k, type(v))
        s = repr(v)
        print(s[:1500])
