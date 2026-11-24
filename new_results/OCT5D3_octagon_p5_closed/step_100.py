import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from fractions import Fraction as F
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod, k32_value
from q2_lll import fit

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
Kred = pickle.load(open("kernel_red.pkl","rb"))['Kred']
mixed = pickle.load(open("mixed61_vals.pkl","rb"))
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
    atoms = list(zip(names, allvals))
    return R, atoms
R1100, atoms1100 = setup(1100)
c = fit(R1100, atoms1100, cexp=int(1100*0.92), nbest=1)
rel1100 = c[0][2]
print("1100 relation height %s residual %s" % (mp.nstr(mp.mpf(c[0][0]),3), mp.nstr(c[0][1],3)))
pickle.dump({k: (str(v[0]), str(v[1])) for k, v in rel1100.items()}, open("rel1100_frac.pkl","wb"))
print("saved rel1100_frac.pkl; entries:", len(rel1100))
