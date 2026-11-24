import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_atoms import W3_BASIS, W4_BASIS, prod
from q2_lll import fit

def load(dps):
    mp.mp.dps = dps
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    Wt = pickle.load(open("wparts_true_1650.pkl","rb"))
    return A, Wt

def run(level, basis, tag, precisions=(900, 1500)):
    outs = []
    for D in precisions:
        A, Wt = load(D)
        T = mp.mpf(Wt[f'Wt{level}'])
        atoms = [(nm, prod(A, sp)) for nm, sp in basis]
        c = fit(T, atoms, cexp=int(D*0.92), nbest=2)
        outs.append(c)
        print(f"--- {tag} at {D} ---")
        for h, r, rel, v0 in c:
            print("  height %s residual %s" % (mp.nstr(mp.mpf(h),3), mp.nstr(r,3)))
    r0, r1 = outs[0][0], outs[1][0]
    same = r0[2] == r1[2]
    print(f">>> {tag}: identical at both: {same}")
    if same:
        print("    REL:", {k:(str(v[0]),str(v[1])) for k,v in r0[2].items()})
    print()
    return outs

print("W3_BASIS:", len(W3_BASIS), "atoms;  W4_BASIS:", len(W4_BASIS), "atoms")
o3 = run(3, W3_BASIS, "Wt3")
