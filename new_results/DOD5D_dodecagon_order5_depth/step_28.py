import pickle, sympy as sp
D = pickle.load(open("<project>/research/08_theory/E5_edge_law/steps/step_02_full.pkl","rb"))
print(D.keys())
E5s = D['E5']
print(type(E5s), len(E5s))
print(E5s[:1500])
