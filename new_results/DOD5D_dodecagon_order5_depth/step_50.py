import pickle
bs = pickle.load(open("<project>/research/08_theory/OCT5D2_octagon_depth_test/block_symbols.pkl","rb"))
print(type(bs), bs.keys() if isinstance(bs, dict) else len(bs))
for k, v in bs.items():
    print(k, type(v), len(v) if hasattr(v, "__len__") else "")
    if isinstance(v, dict):
        items = list(v.items())[:3]
        for kk, vv in items: print("   ", kk, type(vv), vv)
