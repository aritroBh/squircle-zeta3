import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod
from q2_lll import fit

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
kr = pickle.load(open("kernel_red.pkl","rb"))
Kred = kr['Kred']
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
vis = [(lab, spec_map[lab]) for lab in labels]

def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*mp.sqrt(2)

def atom_vals(D):
    mp.mp.dps = D
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    V = [prod(A, spec) for lab, spec in vis]
    return A, V

def compute_R(A, V, D):
    mp.mp.dps = D
    Wt5 = mp.mpf(pickle.load(open("wparts_true_1650.pkl","rb"))['Wt5'])
    classical = sum(q2f(x0[j])*V[j] for j in range(len(V)))
    Cval = sum(q2f(cf)*prod(A, spec) for lab, spec, cf in C_SPECS)
    return Wt5 - classical - Cval

def kernel_combo_values(V, D):
    mp.mp.dps = D
    r2 = mp.sqrt(2)
    out = []
    for row in Kred:
        v = mp.mpf(0)
        for j in range(len(labels)):
            p, q = row[2*j], row[2*j+1]
            if p or q: v += (mp.mpf(p) + mp.mpf(q)*r2)*V[j]
        out.append(v)
    return out

t0 = time.time()
outs = []
for D in (800, 1500):
    A, V = atom_vals(D)
    R = compute_R(A, V, D)
    kvals = kernel_combo_values(V, D)
    atoms = [(nm, prod(A, sp)) for nm, sp in INVISIBLE5] + [(f"K{i+1}", v) for i, v in enumerate(kvals)]
    cands = fit(R, atoms, cexp=int(D*0.92), nbest=3)
    outs.append(cands)
    print(f"--- augmented closure at {D} digits (%.0fs) ---" % (time.time()-t0))
    for h, r, rel, v0 in cands[:2]:
        print("  height %s residual %s" % (mp.nstr(mp.mpf(h),3), mp.nstr(r,3)))
r800, r1500 = outs[0][0], outs[1][0]
print("identical at both:", r800[2] == r1500[2])
if r800[2] == r1500[2]:
    print("REL:", {k:(str(v[0]),str(v[1])) for k,v in r800[2].items()})
print("total %.1fs" % (time.time()-t0))
