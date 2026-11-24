import pickle, sympy as sp
base = '<project>/research/08_theory/'
Jloc = pickle.load(open(base+'OCT6_octagon_closed/oct5/Jloc.pkl','rb'))
print("ALL Jloc entries:")
for k in sorted(Jloc.keys()):
    print(" ", k, "=", Jloc[k], "  syms:", Jloc[k].free_symbols if hasattr(Jloc[k],'free_symbols') else 'n/a')
