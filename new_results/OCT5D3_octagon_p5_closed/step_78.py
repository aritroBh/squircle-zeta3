import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
import mpmath as mp
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod
from q2_lll import fit
from qsqrt2 import QS2 as Q

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
specs = atom_specs()
spec_map = {lab: (kind, spec) for lab, kind, spec in specs}
vis = [(lab, spec_map[lab][1]) for lab in labels]

def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*mp.sqrt(2)

def compute_R(D):
    mp.mp.dps = D
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    Wt5 = mp.mpf(pickle.load(open("wparts_true_1650.pkl","rb"))['Wt5'])
    classical = mp.mpf(0)
    for lab, spec in vis:
        classical += q2f(x0[labels.index(lab)])*prod(A, spec)
    Cval = mp.mpf(0)
    for lab, spec, cf in C_SPECS:
        Cval += q2f(cf)*prod(A, spec)
    R = Wt5 - classical - Cval
    return A, R

t0=time.time()
outs=[]
for D in (800, 1500):
    A, R = compute_R(D)
    atoms = [(nm, prod(A, sp)) for nm, sp in INVISIBLE5]
    cands = fit(R, atoms, cexp=int(D*0.92), nbest=4)
    outs.append((R, cands))
    print(f"--- closure R at {D} digits ---  R = {mp.nstr(R, 12)}")
    for h, r, rel, v0 in cands[:3]:
        print("  height %s residual %s" % (mp.nstr(mp.mpf(h),3), mp.nstr(r,3)))
print("total %.1fs" % (time.time()-t0))
r800, r1500 = outs[0][1][0], outs[1][1][0]
print("identical at both:", r800[2] == r1500[2])
if r800[2] == r1500[2]:
    print("REL:", {k:(str(v[0]),str(v[1])) for k,v in r800[2].items()})
