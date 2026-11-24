import pickle
Jtab = pickle.load(open("<project>/research/08_theory/K1_nonlocal_law/Jtab_t.pkl",'rb'))
E3g = pickle.load(open("<project>/research/08_theory/K1_nonlocal_law/E3_general.pkl",'rb'))
print("Jtab keys:", sorted(Jtab.keys()))
print("n keys:", len(Jtab))
print("E3g keys:", E3g.keys())
print("nl count:", len(E3g['nl']), "loc count:", len(E3g['loc']))
print("sample nl:", E3g['nl'][:3])
