import pickle, os
base = "<project>/research/08_theory/"
D1, D2 = base+"DOD1_dodecagon/", base+"DOD2_dodecagon_finish/"
print(sorted(os.listdir(D1)))
ac = pickle.load(open(D2+"atom_coefs_pq.pkl","rb"))
print("\natom_coefs keys sample:", list(ac.items())[:4], "… total", len(ac))
xv = pickle.load(open(D2+"Xi_verified_620.pkl","rb"))
print("\nXi_verified type:", type(xv), xv if not hasattr(xv,'keys') else xv.keys())
g = pickle.load(open(D2+"gdswaps620.pkl","rb"))
print("\ngdswaps keys:", g.keys() if hasattr(g,'keys') else type(g))
for k in g: 
    print(" ", k, "n=", len(g[k]), "sample:", list(g[k].items())[:1])
