import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod, k32_value
from q2_lll import fit

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
Kred = pickle.load(open("kernel_red.pkl","rb"))['Kred']
mixed = pickle.load(open("mixed61_vals.pkl","rb"))
keep = pickle.load(open("basis_keep.pkl","rb"))['keep']
names = pickle.load(open("basis_keep.pkl","rb"))['names']
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
    allvals = ([prod(A, sp) for nm, sp in INVISIBLE5] + [k32_value()]
               + [mp.mpf(val) for val, desc in mixed] + kvals)
    atoms = [(names[i], allvals[i]) for i in keep]
    return R, atoms

outs = []
for D in (1100, 1500):
    R, atoms = setup(D)
    c = fit(R, atoms, cexp=int(D*0.92), nbest=2)
    outs.append(c)
    print(f"--- reduced 39-atom basis at {D} ---")
    for h, r, rel, v0 in c:
        print("  height %s residual %s" % (mp.nstr(mp.mpf(h),3), mp.nstr(r,3)))
r1, r2_ = outs[0][0], outs[1][0]
print("identical at both:", r1[2] == r2_[2])
if r1[2] == r2_[2]:
    print("ACCEPTED relation:")
    for k, val in r1[2].items():
        print("   %-18s %s" % (k, val))
    pickle.dump({'rel': r1[2], 'height': str(r1[0]), 'names': names, 'keep': keep}, open("closure_accepted.pkl","wb"))
