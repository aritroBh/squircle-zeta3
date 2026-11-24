import sys, importlib, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))

def hip_str(H):
    # +i principal part of higher part H: H[1] = [e1..em] coeff of (v-i)^{-j}
    if not H[1]: return "0"
    return "[" + ", ".join(f"{e.sympy()}" for e in H[1]) + "]   (orders 1..%d)"%len(H[1])

for key in [(2,0),(1,1),(3,0),(2,1),(4,0),(3,1),(2,2)]:
    print(f"===== H(W_{key})  +i principal parts by level k")
    for k in range(sum(key), 11):
        rec = record.get((k,key))
        if rec is None: continue
        print(f"  k={k}: {hip_str(rec['H'])}")
