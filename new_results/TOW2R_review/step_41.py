import pickle, os
p = '<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl'
data = pickle.load(open(p,'rb'))
print(type(data))
if isinstance(data, dict):
    ks = list(data.keys())[:10]
    print('keys sample:', ks)
    k0 = ks[0]
    print('value type:', type(data[k0]))
    print(str(data[k0])[:500])
