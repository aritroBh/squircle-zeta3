import pickle, mpmath as mp, os
base = "<project>/research/08_theory/"
D3 = base+"DOD3_dodecagon_symbol/"
mp.mp.dps = 620
rv = pickle.load(open(D3+"R_all12_vals620.pkl","rb"))
print("triple: |Rv - Rq|")
for k in sorted(rv['Rv']):
    d = abs(rv['Rv'][k]-rv['Rq'][k])
    print("  ", k, mp.nstr(d,3))
# find DOD2 step that verifies T3a_cf vs quadrature
D2 = base+"DOD2_dodecagon_finish/"
for f in sorted(os.listdir(D2)):
    if f.startswith("step_") and f.endswith(".py"):
        txt = open(D2+f).read()
        if "5.232e-150" in txt or ("T3a_cf" in txt and "quad" in txt):
            print("=====", f)
            print(txt[:4000]); break
