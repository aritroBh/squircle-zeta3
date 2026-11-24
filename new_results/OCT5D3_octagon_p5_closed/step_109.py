import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_atoms import W4_BASIS, prod
from q2_lll import fit

def load(dps):
    mp.mp.dps = dps
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    Wt = pickle.load(open("wparts_true_1650.pkl","rb"))
    return A, Wt

t0 = time.time()
outs = []
for D in (900, 1500):
    A, Wt = load(D)
    T = mp.mpf(Wt['Wt4'])
    atoms = [(nm, prod(A, sp)) for nm, sp in W4_BASIS]
    c = fit(T, atoms, cexp=int(D*0.92), nbest=2)
    outs.append(c)
    print(f"--- Wt4 at {D} ({time.time()-t0:.0f}s) ---")
    for h, r, rel, v0 in c:
        print("  height %s residual %s" % (mp.nstr(mp.mpf(h),3), mp.nstr(r,3)))
r0, r1 = outs[0][0], outs[1][0]
same = r0[2] == r1[2]
print(">>> Wt4: identical at both:", same)
if same:
    print("    REL:", {k:(str(v[0]),str(v[1])) for k,v in r0[2].items()})
print("total %.0fs" % (time.time()-t0))
