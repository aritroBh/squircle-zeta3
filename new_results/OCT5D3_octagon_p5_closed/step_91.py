import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod, k32_value
from q2_lll import fit
try:
    import fpylll
    print("fpylll available:", fpylll.__version__ if hasattr(fpylll,'__version__') else 'yes')
except ImportError:
    print("fpylll NOT available, using flint fmpz_mat.lll")

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
Kred = pickle.load(open("kernel_red.pkl","rb"))['Kred']
mixed = pickle.load(open("mixed61_vals.pkl","rb"))   # [(val_str, desc)]
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
vis = [(lab, spec_map[lab]) for lab in labels]

def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*mp.sqrt(2)

def setup(D):
    mp.mp.dps = D
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    V = [prod(A, spec) for lab, spec in vis]
    Wt5 = mp.mpf(pickle.load(open("wparts_true_1650.pkl","rb"))['Wt5'])
    classical = sum(q2f(x0[j])*V[j] for j in range(len(V)))
    Cval = sum(q2f(cf)*prod(A, spec) for lab, spec, cf in C_SPECS)
    R = Wt5 - classical - Cval
    r2 = mp.sqrt(2)
    kvals = []
    for row in Kred:
        v = mp.mpf(0)
        for j in range(len(labels)):
            p, q = row[2*j], row[2*j+1]
            if p or q: v += (mp.mpf(p) + mp.mpf(q)*r2)*V[j]
        kvals.append(v)
    atoms = ([(nm, prod(A, sp)) for nm, sp in INVISIBLE5]
             + [("K32", k32_value())]
             + [(f"MIX{i+1}", mp.mpf(val)) for i, (val, desc) in enumerate(mixed)]
             + [(f"K{i+1}", v) for i, v in enumerate(kvals)])
    return R, atoms

t0 = time.time()
outs = []
for D in (900, 1500):
    R, atoms = setup(D)
    cands = fit(R, atoms, cexp=int(D*0.92), nbest=3)
    outs.append((R, cands))
    print(f"--- full augmented closure at {D} digits ({time.time()-t0:.0f}s), {len(atoms)} atoms ---")
    for h, r, rel, v0 in cands[:2]:
        print("  height %s residual %s" % (mp.nstr(mp.mpf(h),3), mp.nstr(r,3)))
r900, r1500 = outs[0][1][0], outs[1][1][0]
print("identical at both:", r900[2] == r1500[2])
if r900[2] == r1500[2]:
    print("REL:", {k:(str(v[0]),str(v[1])) for k,v in r900[2].items()})
print("total %.1fs" % (time.time()-t0))
