import pickle
Jtab = pickle.load(open('research/08_theory/K1d_engine_fixed/Jtab_t.pkl','rb'))
print(type(Jtab), len(Jtab))
ks = sorted(Jtab.keys())
print(ks[:12])
import sympy as sp
print(Jtab[ks[0]])
Jloc = pickle.load(open('research/08_theory/OCT6_octagon_closed/oct5/Jloc.pkl','rb'))
print(type(Jloc), len(Jloc))
ks2 = sorted(Jloc.keys())[:12]
print(ks2)
print(Jloc[ks2[0]])