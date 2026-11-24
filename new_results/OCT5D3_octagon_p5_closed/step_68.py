import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 800   # set BEFORE creating Afull
from oct5d3_atoms import W1_BASIS, W2_BASIS, W3_BASIS, prod
from q2_lll import fit
Afull = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
Wt = pickle.load(open("wparts_true_1650.pkl","rb"))
T = mp.mpf(Wt['Wt1'])
atoms = [(nm, prod(Afull, sp)) for nm, sp in W1_BASIS]
cands = fit(T, atoms, cexp=735, nbest=3)
for h, r, rel, v0 in cands:
    print("height %s residual %s" % (mp.nstr(mp.mpf(h),3), mp.nstr(r,3)))
    print("   ", {k:(str(v[0]),str(v[1])) for k,v in rel.items()})
