import pickle, os
base = '<project>/research/08_theory'
work = f'{base}/TOW3_tower_depth6'
H = pickle.load(open(f'{base}/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
meta = H.get('_meta'); print('meta', meta)
keys = [k for k in H if k != '_meta']
degs = sorted(set(p+q for p,q in keys))
print('num keys', len(keys), 'degrees', degs[0], '..', degs[-1])
# check a sample key structure
k0 = (5,2)
print('key', k0, 'deltas', sorted(H[k0].keys()))
print('delta 5 levels:', sorted(H[k0][5].keys())[:10], '...', sorted(H[k0][5].keys())[-3:])
print('sample val', H[k0][5][7])
# max degree with delta=6 available and k=d available
ok = []
for key in keys:
    d = sum(key)
    for delta in [5,6]:
        if delta in H[key] and d in H[key][delta] and 2*d-4 >= delta:
            ok.append((key,delta))
print('cells with rho(5):', sum(1 for _,dl in ok if dl==5), ' rho(6):', sum(1 for _,dl in ok if dl==6))
