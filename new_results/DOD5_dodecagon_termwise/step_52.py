import pickle
base = "<project>/research/08_theory/"
lll = pickle.load(open(base+"DOD1_dodecagon/lll102_rels.pkl","rb"))
print(lll.keys())
names = lll['names']
cats = lll['cats']
print({k: len(v) for k,v in cats.items()})
# show T3p relations
for (vc, resid, nrm) in cats.get('T3p', [])[:6]:
    sup = [(names[q], vc[q]) for q in range(len(names)) if vc[q]!=0]
    print("T3p-rel maxc=%s resid=%s: %s"%(nrm, resid[:12], sup[:14]))
