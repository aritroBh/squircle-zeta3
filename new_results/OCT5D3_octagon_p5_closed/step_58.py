import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_atoms import W1_BASIS, W2_BASIS, W3_BASIS, prod
from q2_lll import fit

Afull = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
Wt = pickle.load(open("wparts_true_1650.pkl","rb"))

def run(level, basis, tag, precisions=(800, 1500), extra=()):
    outs = []
    for D in precisions:
        mp.mp.dps = D
        T = mp.mpf(Wt[f'Wt{level}'])
        atoms = [(nm, prod(Afull, sp)) for nm, sp in basis] + list(extra)
        cands = fit(T, atoms, cexp=int(D*0.92), nbest=3)
        outs.append(cands)
        print(f"--- {tag} at {D} digits ---")
        for h, r, rel, v0 in cands[:2]:
            print("  height %s residual %s" % (mp.nstr(mp.mpf(h), 3), mp.nstr(r, 3)))
    r0, r1 = outs[0][0], outs[1][0]
    same = r0[2] == r1[2]
    print(f">>> {tag}: identical at both: {same}; residuals {mp.nstr(r0[1],2)}, {mp.nstr(r1[1],2)}")
    if same:
        print("    REL:", {k: (str(v[0]), str(v[1])) for k, v in r0[2].items()})
    print()
    return outs

mp.mp.dps = 100
extra = [('pi', mp.pi)]
o1 = run(1, W1_BASIS, "Wt1 (+pi)", extra=extra)
