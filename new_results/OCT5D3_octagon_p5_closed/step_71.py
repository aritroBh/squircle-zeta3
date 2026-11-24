import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_atoms import W3_BASIS, prod
from q2_lll import fit

def load(dps):
    mp.mp.dps = dps
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    Wt = pickle.load(open("wparts_true_1650.pkl","rb"))
    return A, Wt
for D in (800, 1500):
    A, Wt = load(D)
    T = mp.mpf(Wt['Wt3'])
    atoms = [(nm, prod(A, sp)) for nm, sp in W3_BASIS]
    cands = fit(T, atoms, cexp=int(D*0.92), nbest=1)
    h, r, rel, v0 = cands[0]
    print(D, {k: (str(v[0]), str(v[1])) for k, v in sorted(rel.items())})
