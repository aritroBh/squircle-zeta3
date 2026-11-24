import pickle
base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
sf = pickle.load(open(base+"symbol_fit2_solution.pkl","rb"))
print(type(sf))
if isinstance(sf, dict):
    for k,v in sf.items():
        print("KEY:", k, type(v), str(v)[:300])
elif isinstance(sf, (list,tuple)):
    print("len", len(sf))
    for i,v in enumerate(sf[:10]):
        print(i, type(v), str(v)[:300])
