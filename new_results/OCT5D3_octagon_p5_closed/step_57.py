import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_atoms import W1_BASIS, W2_BASIS, W3_BASIS, prod
from q2_lll import fit
from fractions import Fraction as F

Afull = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
Wt = pickle.load(open("wparts_true_1650.pkl","rb"))

def run(level, basis, tag, precisions=(800, 1500)):
    results = []
    for D in precisions:
        mp.mp.dps = D
        T = mp.mpf(Wt[f'Wt{level}'])
        atoms = [(nm, prod(Afull, sp)) for nm, sp in basis]
        cands = fit(T, atoms, cexp=int(D*0.92), nbest=3)
        results.append(cands)
        print(f"--- {tag} at {D} digits ---")
        for h, r, rel, v0 in cands[:2]:
            print("  height %.2e residual %s" % (h, mp.nstr(r, 3)))
    # compare best relations
    return results

for lvl, basis, tag in [(1, W1_BASIS, "Wt1"), (2, W2_BASIS, "Wt2"), (3, W3_BASIS, "Wt3")]:
    res = run(lvl, basis, tag)
    r800, r1500 = res[0][0], res[1][0]
    same = r800[2] == r1500[2]
    print(f">>> {tag}: best relations identical at 800/1500: {same}; heights {r800[0]:.2e}, {r1500[0]:.2e}; residuals {mp.nstr(r800[1],2)}, {mp.nstr(r1500[1],2)}")
    if same:
        print("    REL:", {k: (str(v[0]), str(v[1])) for k, v in r800[2].items()})
    print()
