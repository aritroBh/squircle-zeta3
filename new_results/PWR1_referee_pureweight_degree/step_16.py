import pickle
e_rm = pickle.load(open('../K3_fourier_all_orders/e_rm.pkl','rb'))
print(type(e_rm))
if isinstance(e_rm, dict):
    ks = sorted(e_rm.keys())[:40]
    for k in ks: print(k, e_rm[k])
