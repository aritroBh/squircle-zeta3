import pickle
W='<project>/research/08_theory/WTS_weight_k_structure/theory'
for k in range(6,14):
    K=pickle.load(open(f'{W}/kern_{k}.pkl','rb'))
    print(k,[ (key,K[key]['const']) for key in [(6,0),(5,1),(4,2),(3,3)] if any(x!=0 for x in K[key]['const'])] or 'all const zero', 'ipos/ineg empty:', all(not K[key]['ipos'] and not K[key]['ineg'] for key in [(6,0),(5,1),(4,2),(3,3)]))
