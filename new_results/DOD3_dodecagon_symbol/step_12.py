import os, pickle
base = os.path.dirname(os.getcwd())
d2 = os.path.join(base,'DOD2_dodecagon_finish')
with open(os.path.join(d2,'atom_coefs_pq.pkl'),'rb') as f:
    ac = pickle.load(f)
print("atom_coefs keys:")
for k,v in ac.items():
    print('  ', k, v)
