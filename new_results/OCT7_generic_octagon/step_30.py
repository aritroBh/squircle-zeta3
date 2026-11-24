import pickle, sympy as sp
base="research/08_theory"
Jloc = pickle.load(open(base+"/OCT5_octagon_exact/Jloc.pkl","rb"))
print("Jloc keys (sample):", list(Jloc.keys())[:8], "... total", len(Jloc))
k = list(Jloc.keys())[0]
print("sample:", k, "=", Jloc[k])
Jtab = pickle.load(open(base+"/K1d_engine_fixed/Jtab_t.pkl","rb"))
print("\nJtab keys count:", len(Jtab), "sample keys:", list(Jtab.keys())[:6])
k2 = list(Jtab.keys())[0]
print("sample:", k2, "=", Jtab[k2])
import sympy as sp
print("\nfree symbols in Jtab sample:", sorted([str(s) for s in Jtab[k2].free_symbols]))
