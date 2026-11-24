import pickle, sympy as sp
base = '<project>/research/08_theory/'
Jloc = pickle.load(open(base+'OCT6_octagon_closed/oct5/Jloc.pkl','rb'))
Jtab = pickle.load(open(base+'K1d_engine_fixed/Jtab_t.pkl','rb'))
print("Jloc keys:", sorted(Jloc.keys())[:20])
print("Jloc sample:", {k: Jloc[k] for k in sorted(Jloc.keys())[:4]})
print()
print("Jtab keys count:", len(Jtab), " sample keys:", sorted(Jtab.keys())[:8])
# symbols used in Jtab
allsyms=set()
for k,v in Jtab.items():
    allsyms |= v.free_symbols
print("Jtab symbols:", allsyms)
# are coefficients rational? check a couple
for k in sorted(Jtab.keys())[:3]:
    print("Jtab",k,"=",Jtab[k])
print()
allj=set()
for k,v in Jloc.items(): allj|=v.free_symbols
print("Jloc symbols:", allj)
