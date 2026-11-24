
import pickle, os
K3 = '<project>/research/08_theory/K3_fourier_all_orders/'
for f in ['e_rm.pkl','e_rm_ext.pkl','Qpolys.pkl','Dtilde_relations.pkl']:
    with open(K3+f,'rb') as fh:
        d = pickle.load(fh)
    print('='*20, f, type(d))
    print(d)
