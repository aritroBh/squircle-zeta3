import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 1650
from oct5d3_basis import atom_specs
from oct5d3_atoms import prod, li32_xy
t0 = time.time()
kd = pickle.load(open("kernel61.pkl","rb"))
K61, labels, ncl = kd['K61'], kd['labels'], kd['nclassical']
A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
V = [prod(A, spec_map[lab]) for lab in labels[:ncl]]
r2 = mp.sqrt(2); b = r2-1
LI32V = {
 'Li32(b,1)': li32_xy(b, 1), 'Li32(1/sqrt2,1)': li32_xy(1/r2, 1),
 'Li32(-b,1)': li32_xy(-b, 1), 'Li32(b,-1)': li32_xy(b, -1),
 'Li32(-b,-1)': li32_xy(-b, -1), 'Li32(1/sqrt2,-1)': li32_xy(1/r2, -1),
 'Li32(1-1/sqrt2,1)': li32_xy(1-1/r2, 1), 'Li32(2-sqrt2,1)': li32_xy(2-r2, 1),
 'Li32(-b^2,1)': li32_xy(-b*b, 1),
}
def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
mixed = []
for v in K61:
    if any(v[ncl+j] != 0 for j in range(9)):
        val = mp.mpf(0)
        for j in range(ncl):
            if v[j] != 0: val += q2f(v[j])*V[j]
        for j in range(9):
            if v[ncl+j] != 0: val += q2f(v[ncl+j])*LI32V[labels[ncl+j]]
        desc = [(labels[j], str(v[j])) for j in range(len(labels)) if v[j] != 0]
        mixed.append((val, desc))
print("mixed kernel constants (%.1fs):" % (time.time()-t0))
for val, desc in mixed:
    print("  %s   from %s" % (mp.nstr(val, 25), desc[:4]))
pickle.dump([(mp.nstr(val, 1620), desc) for val, desc in mixed], open("mixed61_vals.pkl","wb"))
print("saved", len(mixed))
