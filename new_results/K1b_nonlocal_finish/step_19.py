import pickle, sympy as sp
E3g = pickle.load(open("<project>/research/08_theory/K1_nonlocal_law/E3_general.pkl",'rb'))
allsyms = set()
for key, co in E3g['nl'] + E3g['loc']:
    allsyms |= co.free_symbols
print("all symbols in coefficients:", sorted(str(s) for s in allsyms))
# print all 40 nonlocal monomials with coefficients
for (n,a,b,c), co in sorted(E3g['nl']):
    print((n,a,b,c), ":", co)
