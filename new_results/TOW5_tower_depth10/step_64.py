import pickle
work = '<project>/research/08_theory/TOW5_tower_depth10'
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6, 14)}
keys6 = [(6,0),(5,1),(4,2),(3,3)]
for key in keys6:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    print(key, "degree-6 kernel real poles:", sorted(poles))
# Also count components per key
tot=0
for key in keys6:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    tot += len(poles)*4
    print(key, len(poles), "poles ->", len(poles)*4, "components")
print("total components:", tot)
